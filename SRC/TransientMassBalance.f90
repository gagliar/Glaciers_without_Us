SUBROUTINE TransientMassBalance( Model,Solver,dt,TransientSimulation )
  
  USE DefUtils
  USE SolverUtils
  USE ElementUtils

IMPLICIT NONE
TYPE(Model_t) :: Model
TYPE(Variable_t), POINTER :: Accumulation,Rad_fact_var,SurfGrad1,Surfgrad2,Time
TYPE(Variable_t), POINTER :: MB,Dens,Firn,SE,Depth,Melting,Refreeze,Raining,PotRad
TYPE(Solver_t), POINTER :: Solver
TYPE(Element_t),POINTER :: Element
TYPE(ValueList_t), POINTER :: SolverParams
INTEGER, POINTER :: NodeIndexes(:)

INTEGER :: n,i,j,k,cont,nb_surf,nb_vert,io,nb_year,nb_day,it,day,YearDay,nb_jour,first_day,nb_line
REAL(KIND=dp) :: f, z, deg_pos,accu_ref,a,accu,melt_local,accu_ice, temp_10m,rain,t_simu,deg_jour,T,dt
REAL(KIND=dp) :: z_precip,grad_accu,grad,z_temp,seuil_precip,Rad_fact_snow,Rad_fact_ice,Rad_fact
REAL(KIND=dp) :: Pfact,temp_correc,surimposed_ice_fact,firn_param,deg1,deg2,precip_correc,sigma,melt

REAL(KIND=dp) :: Sr,rho_w,rho_ice,L_heat,rho_surf,T0,g1,g2,g3,Mean_Temp_Air,x_output,y_output

REAL(KIND=dp) :: fx,fy,slop,asp,S0,dr,lat,L,term1,term2,term3,tau_r,tau_d,tau_b,srad,sinAlpha,R,M,Is
REAL(KIND=dp) :: Ir,Iday,I0,hsr,hs,cos_i,dS,Idiff,reflec,Norm,Year0,year,x,y,yearstartdata,precip_fact,MinAltFact,MaxAltFact

real,dimension(:,:),allocatable :: Mask_l
real,dimension(190,223,3) :: Mask

real (KIND=dp), dimension(:),allocatable :: TempAirMin,TempAirMax,TempAirMoy,Precip,PotRadNodes,FirnNodes,MaskAccu
real (KIND=dp),dimension(365) :: TempAirMean, PrecipMean,TempAirMeanTry,TimeStep

character(LEN=MAX_NAME_LEN) :: filenameMax,filenameMin,filename2,filenameMoy

logical :: first_time=.true.,TransientSimulation, GotIt, PrecipData, node_output=.false.,GotNode=.false.,Output1D,SigmaOK
logical :: OutputFirn,OutputDens,OutputMB,OutputMelting,OutputAccumulation,Execute_steady=.true.
logical :: OutputRefreeze,OutputRad_Fact_var,OutputRaining,OutputPotRad,TransientMB

save first_time,nb_surf,nb_vert,TempAirMin,TempAirMax,TempAirMoy,TempAirMean,nb_day,nb_year,&
&PrecipData,Precip,PrecipMean,PotRadNodes,FirnNodes,TransientMB,Execute_steady,MaskAccu


!===============================================================================
!Initialization=================================================================

if (first_time) then
	first_time=.false.
	
	SolverParams => GetSolverParams()
	
	TransientMB =  GetLogical( SolverParams,'Transient MB', GotIt)
	IF (.NOT. GotIt) TransientMB = .false.

	filenameMin= ListGetString(Model % Constants,'MinAirTemperatureFile', GotIt)
	filenameMax= ListGetString(Model % Constants,'MaxAirTemperatureFile', GotIt)
	filenameMoy= ListGetString(Model % Constants,'MoyAirTemperatureFile', GotIt)
	
	IF (.NOT.GotIt) THEN
      CALL FATAL('Surface Boundary', 'No file for daily air temperature defined in Constants section (AirTemperatureFile=)')
    END IF
	
	filename2= GetString(Model % Constants,'PrecipFile', PrecipData)
	
	IF (.NOT.PrecipData) THEN
		IF (ParEnv % MyPE<1) then
			CALL WARN('Surface Boundary', 'No file for daily precipitation defined in Constants section (PrecipFile=)')
			print*,'Using constant precipitation rate defined by the "Precip" keyword'
		ENDIF
    END IF
	
	open(1,file=filenameMin,status='old')
	nb_day = 0
	do
		read(1,*,iostat=io)
		nb_day = nb_day + 1
		if (io/=0) exit
	enddo
	close(1)
	
	nb_day=nb_day-1
	
	nb_year=floor(nb_day/365.25)
	
	IF (PrecipData) then
		open(1,file=filename2,status='old')
		cont = 0
		do
			read(1,*,iostat=io)
			cont = cont + 1
			if (io/=0) exit
		enddo
		close(1)
		cont=cont-1
	
		IF (cont.ne.nb_day) then
			print*, 'Lenght precip file = ',cont
			print*, 'Lenght temperature file = ', nb_day
			CALL FATAL('Surface MB', 'Precip and temperature files have to be same lenght')
		ENDIF
	
		allocate(Precip(nb_day))
	
		open(1,file=filename2,status='old')
		do i=1,nb_day
			read(1,*) Precip(i)
		enddo
		close(1)
	
	
	
	ENDIF
	
	allocate(TempAirMin(nb_day))
	allocate(TempAirMax(nb_day))
	allocate(TempAirMoy(nb_day))

	open(1,file=filenameMin,status='old')
	do i=1,nb_day
		read(1,*) TempAirMin(i)
	enddo
	close(1)
	
	open(1,file=filenameMax,status='old')
	do i=1,nb_day
		read(1,*) TempAirMax(i)
	enddo
	close(1)
	
	open(1,file=filenameMoy,status='old')
	do i=1,nb_day
		read(1,*) TempAirMoy(i)
	enddo
	close(1)
	
	nb_year=floor(nb_day/365.25)
	

!Get number of surface nodes-----------------------------------
!--------------------------------------------------------------

	Depth => VariableGet( Model % Variables, 'Depth')
	IF ( .not. ASSOCIATED( Depth ) ) THEN
		CALL FATAL('Surface MB','Need variable Depth from Flowdepth Solver')
	ENDIF
	

	cont=0

	DO n=1,model % NumberOfNodes
		if (Depth % Values (Depth % perm (n))==0.0) then
			cont=cont+1
		endif
	ENDDO

	nb_surf=cont
	nb_vert=model % NumberOfNodes/cont
	
	
	allocate(FirnNodes(model % NumberOfNodes))
	allocate(PotRadNodes(model % NumberOfNodes))
	allocate(MaskAccu(model % NumberOfNodes))
	
	
	  
  open(1,file='../Data/MaskAccu.dat')
  
  nb_line = 0
	do
		read(1,*,iostat=io)
		nb_line = nb_line + 1
		if (io/=0) exit
	enddo
	close(1)
	
	nb_line=nb_line-1
	
allocate(Mask_l(nb_line,3))
  
  
  open(1,file='../Data/MaskAccu.dat')
  do i=1,nb_line
    read(1,*) Mask_l(i,1),Mask_l(i,2),Mask_l(i,3)
  enddo
  close(1)

  cont=0
  do i=1,190
    do j=1,223
      cont=cont+1

      Mask(i,j,1)=Mask_l(cont,1)
      Mask(i,j,2)=Mask_l(cont,2)
      Mask(i,j,3)=Mask_l(cont,3)
    enddo
  enddo



DO n=1,model % NumberOfNodes

     x = model % nodes % x (n)
     y = model % nodes % y (n)


k=floor((x-Mask(1,1,1))/40)+1
j=floor((y-Mask(1,1,2))/40)+1

if ((j<=1).or.(j>=190).or.(k<=1).or.(k>=223)) then

MaskAccu(n)=1.0

else


MaskAccu(n)=Mask(j,k,3)*(Mask(j,k+1,1)-x)*(Mask(j+1,k,2)-y)+Mask(j,k+1,3)*(x-Mask(j,k,1))&
&*(Mask(j+1,k,2)-y)+Mask(j+1,k,3)*(Mask(j,k+1,1)-x)*(y-Mask(j,k,2))+Mask(j+1,k+1,3)*&
&(x-Mask(j,k,1))*(y-Mask(j,k,2))
MaskAccu(n)=MaskAccu(n)/40/40

endif

ENDDO
	

	
!end first_time
endif
	
if (Execute_steady) then
	
!===============================================================================
!Get parameter and allocate 1D variables========================================

Depth => VariableGet( Model % Variables, 'Depth')
	IF ( .not. ASSOCIATED( Depth ) ) THEN
		CALL FATAL('Surface MB','Need variable Depth from Flowdepth Solver')
	ENDIF
	
SurfGrad1 => VariableGet( Model % Variables, 'SurfGrad1')
	IF ( .not. ASSOCIATED( SurfGrad1 ) ) THEN
		CALL FATAL('Surface MB','Need variable SurfGrad1 from Flowdepth Solver')
	ENDIF

SurfGrad2 => VariableGet( Model % Variables, 'SurfGrad2')
	IF ( .not. ASSOCIATED( SurfGrad2 ) ) THEN
		CALL FATAL('Surface MB','Need variable SurfGrad2 from Flowdepth Solver')
	ENDIF


Firn => VariableGet( Model % Variables, 'Firn')
	IF ( .not. ASSOCIATED( Firn ) ) THEN
		OutputFirn=.false.
		ELSE
		OutputFirn=.true.
	ENDIF
	
MB => VariableGet( Model % Variables, 'Mass Balance')
	IF ( .not. ASSOCIATED( MB ) ) THEN
		OutputMB=.false.
		ELSE
		OutputMB=.true.
	ENDIF
	
Melting => VariableGet( Model % Variables, 'Melting')
	IF ( .not. ASSOCIATED( Melting ) ) THEN
		OutputMelting=.false.
		ELSE
		OutputMelting=.true.
	ENDIF

Accumulation => VariableGet( Model % Variables, 'Accu')
	IF ( .not. ASSOCIATED( Accumulation ) ) THEN
		OutputAccumulation=.false.
		ELSE
		OutputAccumulation=.true.
	ENDIF
Refreeze => VariableGet( Model % Variables, 'Refreeze')
	IF ( .not. ASSOCIATED( Refreeze ) ) THEN
		OutputRefreeze=.false.
		ELSE
		OutputRefreeze=.true.
	ENDIF
	
Rad_Fact_var => VariableGet( Model % Variables, 'Rad_Fact')
	IF ( .not. ASSOCIATED( Rad_Fact_var ) ) THEN
		OutputRad_Fact_var=.false.
		ELSE
		OutputRad_Fact_var=.true.
	ENDIF
	
Raining => VariableGet( Model % Variables, 'Rain')
	IF ( .not. ASSOCIATED( Raining ) ) THEN
		OutputRaining=.false.
		ELSE
		OutputRaining=.true.
	ENDIF
	
PotRad => VariableGet( Model % Variables, 'PotRad')
	IF ( .not. ASSOCIATED( PotRad ) ) THEN
		OutputPotRad=.false.
		ELSE
		OutputPotRad=.true.
	ENDIF


z_temp=GetConstReal(Model % Constants, "z_temp")
z_precip=GetConstReal(Model % Constants, "z_precip")
seuil_precip=GetConstReal(Model % Constants, "seuil_precip")
surimposed_ice_fact=GetConstReal(Model % Constants, "super_ice")
grad=GetConstReal(Model % Constants, "GradTemp")
Rad_fact_ice=GetConstReal(Model % Constants, "RadFact_ice")
Rad_fact_snow=GetConstReal(Model % Constants, "RadFact_snow")
deg_jour=GetConstReal(Model % Constants, "Deg_jour")
firn_param=GetConstReal(Model % Constants, "firn_param")

grad_accu=GetConstReal(Model % Constants, "GradPrecip")
MaxAltFact=GetConstReal(Model % Constants, "MaxAltFact")
MinAltFact=GetConstReal(Model % Constants, "MinAltFact")

temp_correc=GetConstReal(Model % Constants, "TempCorrec",GotIt)

rho_surf=GetConstReal(Model % Constants, "rho_surf")
rho_w=GetConstReal(Model % Constants, "Water Density")
rho_ice=GetConstReal(Model % Constants, "Ice Density")

	IF (.not.GotIt) THEN
		temp_correc=0.0
	ENDIF


IF (PrecipData) then
	precip_correc=GetConstReal(Model % Constants, "PrecipCorrec",GotIt)
	IF (.not.GotIt) THEN
		precip_correc=1.0
	ENDIF
ELSE
	Pfact=GetConstReal(Model % Constants, "Precip",GotIt)
	IF (.not.GotIt) THEN
		CALL FATAL('Surface MB','No precipition file, need to define mean precipition (Precip = )')
	ENDIF
ENDIF


!==============================================================================
!Compute potential solar radiation=================================================


DO n=1,model % NumberOfNodes
	if (Depth % Values (Depth % perm (n))==0.0) then
	
	
S0 = 1367         
dr= 0.0174532925
lat=GetConstReal(Model % Constants, "Latitude")
reflec=0.6

fx=SurfGrad1 % Values (SurfGrad1 % perm (n))
fy=SurfGrad2 % Values (SurfGrad2 % perm (n))

slop=atan(sqrt(fx**2+fy**2))
asp=atan2(fx,fy)*(-1)
L=lat*dr

term1 = sin(L)*cos(Slop) - cos(L)*sin(Slop)*cos(Asp)
term2 = cos(L)*cos(Slop) + sin(L)*sin(Slop)*cos(Asp)
term3 = sin(Slop)*sin(Asp)

srad=0.0
do i=1,365 
    
I0 = S0 * (1.0 + 0.0344*cos(360.0*dr*(194.0+i)/365.0))
dS = 23.45 * dr* sin(360.0*dr * ( (284.0+i)/365.0 ) )
hsr = acos(-tan(L)*tan(dS))
It=nint(12.0*(1.0+hsr/Pi)-12.0*(1.0-hsr/Pi))
Iday=0
     do j=1,It 
             
        hs=hsr-15.0*dr*j     
        sinAlpha = sin(L)*sin(dS)+cos(L)*cos(dS)*cos(hs)
        M=sqrt(1229.0+((614.0*sinAlpha))**2)-614.0*sinAlpha
        tau_b = 0.56 * (exp(-0.65*M) + exp(-0.095*M))
        tau_d = 0.271-0.294*tau_b
        tau_r = 0.271+0.706*tau_b
        cos_i = (sin(dS)*term1) + (cos(dS)*cos(hs)*term2) + (cos(dS)*term3*sin(hs))
        Is = I0 * tau_b
        R = Is * cos_i
		if (R<0.0) then
			R=0
		endif
        Idiff = I0 * tau_d * cos(Slop)*cos(Slop)/2.0 * sinAlpha
        Ir = I0 * reflec * tau_r * sin(Slop)*sin(Slop)/2.0 * sinAlpha
        R= R + Idiff + Ir
		if (R<0.0) then
			R=0
		endif
         Iday=Iday+R
      enddo
srad = srad + Iday/365.0/24.0
enddo  

PotRadNodes(n)=srad

IF (OutputPotRad) THEN	
PotRad % Values (PotRad % perm (n)) = srad
ENDIF

ENDIF
ENDDO

!Initiate to mean mass balance:

!===============================================================================
!Run Mass Balance model=========================================================

	DO n=1,model % NumberOfNodes
	if (Depth % Values (Depth % perm (n))==0.0) then
	
	
	
		z= model % nodes % z(n)
		accu_ref=0.0
		deg_pos=0.0
		rain=0.0

!===============================================================================
!Get mean steady climate forcing================================================		
		Mean_Temp_Air=0.0
		melt_local=0.0
		
		DO i=1,nb_day
			
			if (PrecipData) then
			Pfact=Precip(i)*365.25*precip_correc
			endif
		
			T=TempAirMoy(i)+grad*(z_temp-z)+temp_correc
			Mean_Temp_Air=Mean_Temp_Air+T/nb_day
			
			precip_fact=max((1.0+(z-z_precip)*grad_accu)*MaskAccu(n),MinAltFact)
			precip_fact=min(precip_fact,MaxAltFact)
    
			melt=(T)*deg_jour+rad_fact_snow*PotRadNodes(n)
			if (melt>0) then
				melt_local=melt_local+melt
			endif
	
			if (T<=seuil_precip) then
				accu_ref=accu_ref+Pfact/365.25*precip_fact
			endif

			if (T>seuil_precip) then
				rain=rain+Pfact/365.25*precip_fact
			endif

		ENDDO


melt_local=melt_local/nb_year
accu_ref=accu_ref/nb_year
rain=rain/nb_year

!===============================================================================
!Compute MB (Gilbert et al., 2016)==============================================

accu=accu_ref
accu_ice = min(surimposed_ice_fact*accu,melt_local)

a=-melt_local-accu_ice*(1.0+1.0/((1.0-rho_surf/rho_ice)*rho_w/rho_surf))+accu

if (a>=0.0) then
	Rad_fact=Rad_fact_snow
else
	Rad_fact=Rad_fact_ice-(Rad_fact_ice-Rad_fact_snow)*&
	&(accu-accu_ice*(1.0+1.0/((1.0-rho_surf/rho_ice)*rho_w/rho_surf)))/melt_local
endif

melt_local=0.0
		
DO i=1,nb_day
	T=TempAirMoy(i)+grad*(z_temp-z)+temp_correc
	melt=(T)*deg_jour+rad_fact*PotRadNodes(n)
	if (melt>0) then
			melt_local=melt_local+melt
	endif
ENDDO

melt_local=melt_local/nb_year



!===============================================================================
!Compute Firn Thickness=========================================================

FirnNodes(n) = a*firn_param !Firn % values (Firn % perm(n)) + a*dt - Firn % values (Firn % perm(n))*dt/firn_param


if (FirnNodes(n)<0.0) then
	FirnNodes(n)= 0.0
endif

if (OutputFirn) then
Firn % values (Firn % perm(n)) = FirnNodes(n)
endif

if (OutputMB) then
if (FirnNodes(n)>=0.0) then
MB % values (MB % perm(n)) = (accu_ice+accu-melt_local)/(rho_surf/rho_w)
else
MB % values (MB % perm(n)) = (accu_ice+accu-melt_local)/(rho_ice/rho_w)
endif
endif

if (OutputMelting) then
Melting % values (Melting % perm(n)) = melt_local
endif
if (OutputRaining) then
Raining % values (Raining % perm(n)) = rain
endif
if (OutputAccumulation) then
Accumulation % values (Accumulation % perm(n)) = accu
endif
if (OutputRad_Fact_var) then
Rad_fact_var % values (Rad_fact_var % perm(n)) = Rad_fact
endif

endif !Surface Node
ENDDO



endif !Execute Steady
	


if (TransientMB) then

Execute_steady=.false.

Depth => VariableGet( Model % Variables, 'Depth')
	IF ( .not. ASSOCIATED( Depth ) ) THEN
		CALL FATAL('Surface MB','Need variable Depth from Flowdepth Solver')
	ENDIF
	
SurfGrad1 => VariableGet( Model % Variables, 'SurfGrad1')
	IF ( .not. ASSOCIATED( SurfGrad1 ) ) THEN
		CALL FATAL('Surface MB','Need variable SurfGrad1 from Flowdepth Solver')
	ENDIF

SurfGrad2 => VariableGet( Model % Variables, 'SurfGrad2')
	IF ( .not. ASSOCIATED( SurfGrad2 ) ) THEN
		CALL FATAL('Surface MB','Need variable SurfGrad2 from Flowdepth Solver')
	ENDIF

Firn => VariableGet( Model % Variables, 'Firn')
	IF ( .not. ASSOCIATED( Firn ) ) THEN
		OutputFirn=.false.
		ELSE
		OutputFirn=.true.
	ENDIF
	
Dens => VariableGet(Model % Mesh % Variables, 'Densi' )
	IF ( .not. ASSOCIATED( Dens ) ) THEN
		OutputDens=.false.
		ELSE
		OutputDens=.true.
	ENDIF

MB => VariableGet( Model % Variables, 'Mass Balance')
	IF ( .not. ASSOCIATED( MB ) ) THEN
		OutputMB=.false.
		ELSE
		OutputMB=.true.
	ENDIF
	
Melting => VariableGet( Model % Variables, 'Melting')
	IF ( .not. ASSOCIATED( Melting ) ) THEN
		OutputMelting=.false.
		ELSE
		OutputMelting=.true.
	ENDIF

Accumulation => VariableGet( Model % Variables, 'Accu')
	IF ( .not. ASSOCIATED( Accumulation ) ) THEN
		OutputAccumulation=.false.
		ELSE
		OutputAccumulation=.true.
	ENDIF
Refreeze => VariableGet( Model % Variables, 'Refreeze')
	IF ( .not. ASSOCIATED( Refreeze ) ) THEN
		OutputRefreeze=.false.
		ELSE
		OutputRefreeze=.true.
	ENDIF
	
Rad_Fact_var => VariableGet( Model % Variables, 'Rad_Fact')
	IF ( .not. ASSOCIATED( Rad_Fact_var ) ) THEN
		OutputRad_Fact_var=.false.
		ELSE
		OutputRad_Fact_var=.true.
	ENDIF
	
Raining => VariableGet( Model % Variables, 'Rain')
	IF ( .not. ASSOCIATED( Raining ) ) THEN
		OutputRaining=.false.
		ELSE
		OutputRaining=.true.
	ENDIF
	
PotRad => VariableGet( Model % Variables, 'PotRad')
	IF ( .not. ASSOCIATED( PotRad ) ) THEN
		OutputPotRad=.false.
		ELSE
		OutputPotRad=.true.
	ENDIF


z_temp=GetConstReal(Model % Constants, "z_temp")
z_precip=GetConstReal(Model % Constants, "z_precip")
seuil_precip=GetConstReal(Model % Constants, "seuil_precip")
surimposed_ice_fact=GetConstReal(Model % Constants, "super_ice")
grad=GetConstReal(Model % Constants, "GradTemp")
Rad_fact_ice=GetConstReal(Model % Constants, "RadFact_ice")
Rad_fact_snow=GetConstReal(Model % Constants, "RadFact_snow")
deg_jour=GetConstReal(Model % Constants, "Deg_jour")
firn_param=GetConstReal(Model % Constants, "firn_param")

grad_accu=GetConstReal(Model % Constants, "GradPrecip")
MaxAltFact=GetConstReal(Model % Constants, "MaxAltFact")
MinAltFact=GetConstReal(Model % Constants, "MinAltFact")

temp_correc=GetConstReal(Model % Constants, "TempCorrec",GotIt)

rho_surf=GetConstReal(Model % Constants, "rho_surf")
rho_w=GetConstReal(Model % Constants, "Water Density")
rho_ice=GetConstReal(Model % Constants, "Ice Density")

	IF (.not.GotIt) THEN
		temp_correc=0.0
	ENDIF

IF (PrecipData) then
	precip_correc=GetConstReal(Model % Constants, "PrecipCorrec",GotIt)
	IF (.not.GotIt) THEN
		precip_correc=1.0
	ENDIF
ELSE
	Pfact=GetConstReal(Model % Constants, "Precip",GotIt)
	IF (.not.GotIt) THEN
		CALL FATAL('Surface MB','No precipition file, need to define mean precipition (Precip = )')
	ENDIF
ENDIF

IF (OutputPotRad) THEN	
PotRad % Values (:) = 0.0
ENDIF
if (OutputMB) then
MB % values (:) = 0.0
endif
if (OutputMelting) then
Melting % values (:) = 0.0
endif
if (OutputRaining) then
Raining % values (:) = 0.0
endif
if (OutputAccumulation) then
Accumulation % values (:) = 0.0
endif
if (OutputRad_Fact_var) then
Rad_fact_var % values (:) = 0.0
endif


year0 = ListGetConstReal( Model % Constants, 'year0')
yearstartdata = ListGetConstReal( Model % Constants, 'yearstartdata')

Time => VariableGet( Model % Variables, 'Time' )
year = Time % Values (1) + year0


nb_jour=NINT(dt*365.25)
nb_jour=MAX(nb_jour,1)

   
first_day = floor((year-yearstartdata)*365.25)
YearDay= floor((year-floor(year))*365.0)+1



!==============================================================================
!Compute potential solar radiation=================================================


do day=first_day,first_day+nb_jour-1


DO n=1,model % NumberOfNodes
	if (Depth % Values (Depth % perm (n))==0.0) then
	
	
S0 = 1367         
dr= 0.0174532925
lat=GetConstReal(Model % Constants, "Latitude")
reflec=0.6

fx=SurfGrad1 % Values (SurfGrad1 % perm (n))
fy=SurfGrad2 % Values (SurfGrad2 % perm (n))

slop=atan(sqrt(fx**2+fy**2))
asp=atan2(fx,fy)*(-1)
L=lat*dr

term1 = sin(L)*cos(Slop) - cos(L)*sin(Slop)*cos(Asp)
term2 = cos(L)*cos(Slop) + sin(L)*sin(Slop)*cos(Asp)
term3 = sin(Slop)*sin(Asp)

srad=0.0
i=YearDay+day-first_day
    
I0 = S0 * (1.0 + 0.0344*cos(360.0*dr*(194.0+i)/365.0))
dS = 23.45 * dr* sin(360.0*dr * ( (284.0+i)/365.0 ) )
hsr = acos(-tan(L)*tan(dS))
It=nint(12.0*(1.0+hsr/Pi)-12.0*(1.0-hsr/Pi))
Iday=0
     do j=1,It 
             
        hs=hsr-15.0*dr*j     
        sinAlpha = sin(L)*sin(dS)+cos(L)*cos(dS)*cos(hs)
        M=sqrt(1229.0+((614.0*sinAlpha))**2)-614.0*sinAlpha
        tau_b = 0.56 * (exp(-0.65*M) + exp(-0.095*M))
        tau_d = 0.271-0.294*tau_b
        tau_r = 0.271+0.706*tau_b
        cos_i = (sin(dS)*term1) + (cos(dS)*cos(hs)*term2) + (cos(dS)*term3*sin(hs))
        Is = I0 * tau_b
        R = Is * cos_i
		if (R<0.0) then
			R=0
		endif
        Idiff = I0 * tau_d * cos(Slop)*cos(Slop)/2.0 * sinAlpha
        Ir = I0 * reflec * tau_r * sin(Slop)*sin(Slop)/2.0 * sinAlpha
        R= R + Idiff + Ir
		if (R<0.0) then
			R=0
		endif
         Iday=Iday+R
      enddo
srad = Iday/24.0
  

PotRadNodes(n)=srad

IF (OutputPotRad) THEN	
do i=1,nb_vert
PotRad % Values (PotRad % perm (n-(i-1)*nb_surf)) = PotRad % Values (PotRad % perm (n-(i-1)*nb_surf)) + srad/nb_jour
enddo
ENDIF
ENDIF
ENDDO




!===============================================================================
!Run Mass Balance model=========================================================

	DO n=1,model % NumberOfNodes
	if (Depth % Values (Depth % perm (n))==0.0) then
	
	
		z= model % nodes % z(n)! + FirnNodes(n)/0.45 - FirnNodes(n)
		rain=0.0
		accu_ref=0.0
		
		

			
			if (PrecipData) then
			Pfact=Precip(day)*365.25*precip_correc
			endif
		
			T=TempAirMoy(day)+grad*(z_temp-z)+temp_correc
			
			precip_fact=max((1.0+(z-z_precip)*grad_accu)*MaskAccu(n),MinAltFact)
			precip_fact=min(precip_fact,MaxAltFact)
									
			if (FirnNodes(n)>0.0) then
				Rad_fact=rad_fact_snow*0.77
				melt=(T)*deg_jour+Rad_fact*PotRadNodes(n)
			else
				Rad_fact=Rad_fact_ice*0.77
				melt=(T)*deg_jour+Rad_fact*PotRadNodes(n)
			endif
			
			if (melt<0) then
				melt=0.0
			endif
	
			if (T<=seuil_precip) then
				accu_ref=Pfact/365.25*precip_fact
			endif

			if (T>seuil_precip) then
				rain=Pfact/365.25*precip_fact
			endif


!===============================================================================
!Compute MB (Gilbert et al., 2016)==============================================

accu=accu_ref


if (OutputMB) then
if (FirnNodes(n)<=0) then
MB % values (MB % perm(n)) = MB % values (MB % perm(n)) + ((accu-melt)/(rho_ice/rho_w))*365.25/nb_jour
else
MB % values (MB % perm(n)) = MB % values (MB % perm(n)) + ((accu-melt)/(rho_surf/rho_w))*365.25/nb_jour
endif
endif
if (OutputMelting) then
Melting % values (Melting % perm(n)) = Melting % values (Melting % perm(n))  + melt*365.25/nb_jour
endif
if (OutputRaining) then
Raining % values (Raining % perm(n)) = Raining % values (Raining % perm(n)) + rain*365.25/nb_jour
endif
if (OutputAccumulation) then
Accumulation % values (Accumulation % perm(n)) = Accumulation % values (Accumulation % perm(n)) + accu*365.25/nb_jour
endif
if (OutputRad_Fact_var) then
Rad_fact_var % values (Rad_fact_var % perm(n)) = Rad_fact_var % values (Rad_fact_var % perm(n)) + Rad_fact/nb_jour
endif



!===============================================================================
!Compute Firn Thickness=========================================================

FirnNodes(n) = FirnNodes(n) + (accu-melt) - FirnNodes(n)/firn_param/365.25

if (FirnNodes(n)<0.0) then
	FirnNodes(n)= 0.0
endif

 if (OutputFirn) then
 do i=1,nb_vert
     Firn % values (Firn % perm(n-(i-1)*nb_surf)) = FirnNodes(n)
 enddo
 endif


endif


ENDDO !Nodes loop
ENDDO !Day loop

ENDIF

 

END SUBROUTINE TransientMassBalance



