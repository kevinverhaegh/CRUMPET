# CRUM reaction class reaction.py
# Separated from CRUM.py by holm10
# Changelog
# 200205 - Separated from CRUM.py #holm10
# 200210 - Updated ADAS extrapolation, tidied up code #holm10
 
class REACTION:


    def __init__(self, name, database, coeffs,typ,S,Tarr=0):
        ''' Creates an reaction object
            __init__(name,database,reactants,fragments,*keys)
    
            name        -   Reaction ID/name (string)
            database    -   Database where reaction is taken from (string)
            reactants   -   List of strings of reactant ID's/handles
            fragments   -   List of strings of fragment ID's/handles
            coeffs      -   Reaction rate coefficients for interpolation/calculation
            typ         -   Reaction type, either of
                            'UE'    -   UEDGE type interpolation
                            'ADAS'  -   ADAS type interpolation
                            'RATE'  -   EIRENE type interpolation, 1D or 2D depending on shape of coeffs
                            'COEFFICIENT'   -   Spontaneous/radiative transition rate
                            'SIGMA' -   Sawada type cross-section fit

            Optional parameters
            Tarr (0)    -   Temperature array required for ADAS interpolation

            
        '''
        from numpy import ones,array
        from scipy.interpolate import interp2d,interp1d
        
        # Store the data required to generate the reaction rates
        self.name=name
        self.coeffs=coeffs
        self.database=database
        self.reactants=S[0]
        self.r_mult=S[2]
        self.fragments=S[1]
        self.f_mult=S[3]
        self.coeffs=array(coeffs)
        self.type=typ
        self.Tarr=array(Tarr)
        for i in range(4):
            if isinstance(S[i],str): S[i]=S[i][4:]

        self.S_r=S[4][0]
        self.S_g=S[4][-2]
        self.S_V=S[4][1]
        self.S_e=S[4][-1]


        # Setup multilplier arrays
        #self.r_mult=ones((len(reactants),))
        #self.f_mult=ones((len(fragments),))
        
        # Separate coefficients from reactant species
        #for i in range(len(self.reactants)):    # Loop through reactants
        #    if '*' in self.reactants[i]:    # If there is a multiplier, separate
        #        self.f_mult[i]=float(reactants[i].split('*')[0])    # Update multiplier
        #        self.reactants[i]=reactants[i].split('*')[1].strip()    # Remove mult from reactant string

        # Separate coefficients from fragment species
        #for i in range(len(self.fragments)):    # Loop through fragments
         #   if '*' in self.fragments[i]:    # If there is a multiplier, separate
         #       self.f_mult[i]=float(fragments[i].split('*')[0])    # Update the multiplier
         #       self.fragments[i]=fragments[i].split('*')[1].strip()    # Remove the multiplier from fragment string


        ''' Make interpolation object if necessary '''
        if self.type=='UE': # UEDGE interpolator
            t=[i for i in range(self.coeffs.shape[0])]
            n=[i for i in range(self.coeffs.shape[1])]

            self.interpolation=interp2d(n,t,self.coeffs,kind='linear') # Create 2D interpolation function

        if self.type=='ADAS': # ADAS interpolator
            self.interpolation=interp1d(self.Tarr,self.coeffs,kind='slinear')   # Create 1D interpolation function

    
    def print_reaction(self):
        ''' Returns a string with the reaction 
            print_reaction()
        '''
        ret='' # Initialize string
        ret+='{}_{}: '.format(self.database,self.name) # Append reaction ID
        for i in range(len(self.reactants)):    # Loop through all reactants
            ret+=(self.r_mult[i]!=1)*'{}*'.format(self.r_mult[i])+'{} '.format(self.reactants[i])+(i+1!=len(self.reactants))*'+ '
        ret+='=> ' # Add reactants to string
        for i in range(len(self.fragments)):    # Loop through all fragmetns
            ret+=(self.f_mult[i]!=1)*'{}*'.format(self.f_mult[i])+'{} '.format(self.fragments[i])+(i+1!=len(self.fragments))*'+ ' # Add fragments to string
        
        return ret 

            
    def rate(self,Te,Ti,E=None,ne=None,omegaj=1):
        ''' Returns the reaction rate at the specified temperature
            rate(Te,Ti,*keys)

            Te          -   Electron temperature for evaluation, used if electron reaction [eV]
            Ti          -   Ion temperature for evaluation, used if proton impact [eV]
            
            Optional parameters
            E (None)    -   Target particle energy, used if proton impact [eV]
            ne (None)   -   Electron density, used for UEDGE rates [cm**3]
            omegaj (1)  -   Statistical weight of ADAS rates 
        '''
        from numpy import log,exp,sqrt,pi,inf,array,log10,logspace,linspace
        from scipy.integrate import quad
        from matplotlib.pyplot import figure

        # Find reactant species
        if 'e' in self.reactants:   T=Te # Electron-mediated reaction, use Te
        elif 'p' in self.reactants: T=Ti # Proton impact reaction, use Ti
        else: T=0   # Photoemission, T not used


        # Extrapolate from minimum to zero



        # Get rate based on self.type
        if self.type=='RATE':
            ''' We have an EIRENE polynomial fit '''
            ret=0
            if T<0.5:   
                Tuse=0.5
                coeff=T/Tuse
            else:       
                Tuse=T
                coeff=1

            if len(self.coeffs.shape)==2:
                ''' T,E fit '''
                for i in range(9):
                    for j in range(9):
                        ret+=self.coeffs[i,j]*(log(Tuse)**i)*(log(E)**j)

            elif len(self.coeffs.shape)==1: 
                ''' T fit '''
                for i in range(9):
                    ret+=self.coeffs[i]*(log(Tuse)**i)
                    
            else:
                print('Unknown fit')

            return coeff*exp(ret)

        elif self.type=='COEFFICIENT':
            ''' Coefficient '''
            return self.coeffs

        elif self.type=='SIGMA':
            ''' SAWADA cross-section '''
            # TODO Extend to general species?
            mH2=2*1.6735575e-27 # Assume H2 is reactant 1
            me=9.10938356e-31   # Assume electron is reactant 2
            ev=1.602e-19        # Helper
            VH2=sqrt((2*E*ev)/mH2)  # Get the H2 velocity
            Ve=sqrt((2*T*ev)/me)    # Get the e- velocity
            mr=(me*mH2)/(me+mH2)    # Reduced mass
            vth=sqrt((2*self.coeffs[0]*ev)/mr)  # Thermal CM speed
            
            def sigma(E,Eth,q0,A,Omega,W,gamma,nu): # Calculate sigma from the Sawada fit based on E
                Psi=(nu!=0)*(1-W/E)**nu+(gamma!=0)*(1-(W/E)**gamma) # Get triplet/singlet Psi
                return (E>=Eth)*q0*(A/W**2)*((W/Eth)**Omega)*Psi    # Perform fit
            
            def R(x,T,coeffs):
                # Integrand function as described in JUEL-3858
                return x*sigma(x*T,*self.coeffs)*exp(-x)

            # Perform integration over velocity space according to JUEL-3858
            return (4/sqrt(pi))*sqrt((T*ev)/(2*me))*quad(R,0,inf,args=(T,self.coeffs))[0]


        elif self.type=='ADAS':
            ''' ADAS fit '''
            if T<self.Tarr[0]:
                Tuse=self.Tarr[0]
                coeff=T/Tuse
            else:
                Tuse=T
                coeff=1
            Tuse=min(Tuse,self.Tarr[-1])
#            T=max(T,self.Tarr[0])
            # TODO: How to deal with extrapolation?
            # TODO: figure out what is implied by the statistical weight omegaj - set =1 for now
            # Return the rate as calculated from the ADAS fit, per the ADAS manual
            return 2.1716e-8*(1/omegaj)*sqrt(13.6048/Tuse)*self.interpolation(Tuse) 
        
        elif self.type=='UE':
            ''' UEDGE fit '''
            # Turn the temperature and density into log-log variables, bounded to the limits
            jt=max(0,min(10*(log10(Te+1e-99)+1.2),60))
            jn=max(0,min(2*(log10(ne)-10),15))
            # Interpolate jt
            
            # If the density is being extrapolated, scale accordingly!
            c=1
            if self.name in ['RECRAD','IONIZRAD']: c=6.242e11
            return self.interpolation(jn,jt)[0]*c

        elif self.type=='APID':

            def sigma(T,coeffs):
                ''' The APID-4 rates are taken from Janev's 1993 IAEA paper, the analytic solutions from Stotler's svlib routine used in DEGAS2 '''
                n=self.coeffs[0]
                I=13.6/n**2

                def expint(k,p):
                    a=[-0.57721566,0.99999193,-0.24991055,0.05519968,-0.00976004,0.00107857]
                    ap=[8.5733287401,18.0590169730,8.6347608925, 0.2677737343 ]
                    b=[9.5733223454,25.6329561486,21.0996530827, 3.9584969228 ]

                    def en(zn,z):
                        return exp(-z)*( 1 + zn/(z+zn)**2 + zn*(zn-2*z)/(z+zn)**4 + zn*(6*z**2-8*zn*z+zn**2)/(z+zn)**6 )/(z+zn)

                    if k==0: 
                        if p<0:
                            return None
                        else:
                            return exp(-p)/p

                    if p<0: 
                        return None
                    elif p==0:
                        if k==1:
                            return None
                        else:
                            return 1/(k-1)

                    elif p<8 or k==1:
                        if p<1:
                            ze1=a[5]
                            for i in [4,3,2,1,0]:
                                ze1=a[i]+ze1*p
                            ze1-=log(p)

                        else:
                            znum=1
                            zden=1
                            for i in range(4):
                                znum=ap[i]+znum*p
                                zden=b[i]+zden*p
                            ze1=exp(-p)*znum/(zden*p)

                        if k==1:
                            return ze1
                        else:
                            zek=ze1
                            for i in range(2,k+1):
                                zek=(exp(-p)-p*zek)/(i-1)
                            return zek
                    
                    else:
                        kp=int(p+0.5)
                        if k<kp and k<=10:
                            if kp>10: kp=10
                            zek=en(kp,p)
                            for i in range(kp-1,k-1,-1): zek=(exp(-p) -i*zek)/p
                            return zek
                    
                        else: return en(k,p)
                    return 'Fell througFell throughh'

                    #return exp(-p)/(p**k)

                # Use set constants to eval sigma
                if n<=3: 
                    A=self.coeffs[1]
                    b=self.coeffs[2:7]
                    zarg=I/T

                    zeint=[]
                    for i in range(1,7): 
                        zeint.append(expint(i,zarg))

                    zmul=[]
                    zmul.append(b[0] + 2*b[1] + 3*b[2] + 4*b[3] + 5*b[4])
                    zmul.append(-2*( b[1] + 3*b[2] + 6*b[3] + 10*b[4] ))
                    zmul.append(3*( b[2] + 4*b[3] + 10*b[4] ))
                    zmul.append(-4*( b[3] + 5*b[4] ))
                    zmul.append(5*b[4])

                    zi1=0
                    for i in range(5): zi1+=zmul[i]*zeint[i+1]


                    zi2=A*zeint[0]
                    zi3=1e-13/(I*T)
                    return 6.692e7*sqrt(T) *zi3*(zi1 + zi2)
      

                # Higher n, evaluate from formulae
                else:

                    zn=n
                    pt=T
                    pin=I
                    zg0 = 0.9935 + 0.2328/zn - 0.1296/zn**2
                    zg1 = -(0.6282 - 0.5598/zn + 0.5299/zn**2) / zn
                    zg2 = -(0.3887 - 1.181/zn + 1.470/zn**2) / zn**2
                    zmul = 32. / (3. * sqrt(3.) * pi)
                    zan = zmul * zn * (zg0/3. + zg1/4. + zg2/5.)  
                    zrn = 1.94 * zn**(-1.57)
                    zb = (4.0 - 18.63/zn + 36.24/zn**2 - 28.09/zn**3) / zn
                    zbn = 2. * zn**2 * (5. + zb) / 3
                    zyn = pin / pt
                    zzn = zrn + zyn
                    ze1y = expint(1,zyn)
                    ze1z = expint(1,zzn)
                    zint1 = zan * (ze1y/zyn - ze1z/zzn) 
                    zxiy = expint(0,zyn) - 2.*ze1y + expint(2,zyn)
                    zxiz = expint(0,zzn) - 2.*ze1z + expint(2,zzn)
                    zint2 = (zbn - zan * log(2.*zn**2)) * (zxiy - zxiz)
                    zint3 = 1.76e-16 * zn**2 * zyn**2
                    ret = 6.692e7 * sqrt(pt) * zint3 * (zint1 + zint2)



                return ret
                    
            def R(x,T,coeffs):
                # Integrand function as described in JUEL-3858
                return x*sigma(x*T,self.coeffs)*exp(-x)

            if self.name.split('_')[0].upper()=='IONIZ':
                ''' APID ionization rate '''

                '''
                f=figure()
                ax=f.add_subplot()
                x=logspace(-2,4,500)
                y=[]
                for i in x:
                    y.append(sigma(i,self.coeffs))
                
                ax.loglog(x,y)
                ax.grid(which='both') 
                ax.set_ylim((1e-10,1e-5))
                '''
                return sigma(T,self.coeffs)

                    

                    
        else:
            print('Unknown type "{}"'.format(self.type))
    



