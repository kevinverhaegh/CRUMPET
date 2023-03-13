#script for calculating effective rates based on CRUMPET

def run_demo(input_crm,iso_mass=1,T_request=1):
    print('Running DEMO with input:')
    print(input_crm)
    import CRUMPET
    import numpy as np

    indx_H2v = np.append(2,np.arange(3,17))

    #initialise CRM
    crm = CRUMPET.Crumpet(input_crm)

    #make Te & ne vectors
    Tev = np.linspace(0.2,10,100)
    Tiv = Tev # Assume Ti=Te
    ne = 1e19*np.ones(np.shape(Tev)) #assume electron density is 1e19 m-3 (should not impact the rates)
    crm.source[2] = 1e-100 #add small source for numerical stability (only needed if reactions that dissociate are included)

    #compute vibrational distribution H2
    fv_H2 = np.zeros([15,len(Tev)])

    for i in range(0,len(Tev)):
        fv_H2[:,i]=crm.steady_state(Tev[i],ne[i],plot=False,dt=True)[indx_H2v]

    #normalise vibrational distribution by dividing the distribution values to the sum of the distribution
    fv_H2 = fv_H2/(np.sum(fv_H2,axis=0)[None,:])

    #Get vibrationally resolved molecular CX rates from H2VIBR

    X=CRUMPET.ratedata.RateData(rates={'H2VIBR' : '/rates/h2vibr.tex', 'AMJUEL' : '/rates/amjuel.tex'})

    #Get rates as function of Te

    def eval_1D(coeff,T):
        o = np.zeros(np.shape(T))
        for i in range(0,len(coeff)):
            o = o + coeff[i]*(np.log(T)**i)
        return 1e-6*np.exp(o)

    vibr_resolved_CX = np.zeros([15,len(Tev)])
    for i in range(0,15):
        vibr_resolved_CX[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'L2'],Tiv/iso_mass)

    #Now use fv_H2 as a weight and sum the total reaction rate to generate the effective rate
    eff_mol_cx = np.sum(vibr_resolved_CX*fv_H2,axis=0)

    #calculate default effective rate Eirene (AMJUEL H.2 3.2.3)
    h3_2_3_tab = eval_1D(X.reactions['AMJUEL']['H.2']['3.2.3'],Tev/iso_mass)

    #calculate effective molecular CX rate using Ichihara rates
    import mat73, os
    ichi_tables = os.path.join(os.path.dirname(__file__), 'MolCX_Ichi.mat')

    ichi_tab = mat73.loadmat(ichi_tables)
    from scipy.interpolate import interp1d
    #make interpolation object
    f_ichi = interp1d(ichi_tab['T'], 1e-6*ichi_tab['MolCX_Ichi_01'], bounds_error=False,
                      fill_value=(1e-6*ichi_tab['MolCX_Ichi_01'][:, 0], 1e-6*ichi_tab['MolCX_Ichi_01'][:, -1])) #interpolate Ichihara tables at EH2=0.1 eV. Nearest neighbour extrapolation

    eff_mol_cx_ichi = np.sum(f_ichi(Tiv/iso_mass)*fv_H2,axis=0)

    #show a plot comparing different effective rates
    import matplotlib.pyplot as plt
    plt.figure()
    plt.loglog(Tev,np.transpose(fv_H2))
    plt.xlabel('Temperature (eV)')
    plt.ylabel('Fractional abundance of vibrational distribution')
    plt.title('Vibr. distribution as function of T')

    plt.figure()
    plt.loglog(Tev,eff_mol_cx,label='Effective rate (NEW)')
    plt.loglog(Tev,h3_2_3_tab,label='Tabulated rate')
    plt.loglog(Tev,eff_mol_cx_ichi,label='Effective rate - Ichihara (NEW)')
    plt.xlabel('Temperature (eV)')
    plt.legend()

    #Compare vibrational distribution to Boltzmann and map to upper state
    indx = np.argmin(np.abs(Tev-T_request))
    plt.figure()
    plt.plot(fv_H2[:,indx]/fv_H2[0,indx],label='Vibr. distribution - GROUND STATE')
    #fit temperature distribution to fv_H2
    import fc_mapping
    from scipy.optimize import curve_fit
    def fit_model(x,a):
        return fc_mapping.vibr_boltzmann(a)[0:len(x)]
    o,_ = curve_fit(fit_model,np.arange(0,15),fv_H2[:,indx]/fv_H2[0,indx],p0=9000,bounds=(100,2e5),method='trf')
    fit_res = fit_model(np.arange(0,15),o)
    plt.plot(fit_res,label='Boltzmann fit T_vib=' + str(o))
    fv_ElecExc_H2 = fc_mapping.get_upper(ground=fv_H2[:,indx])
    plt.plot(fv_ElecExc_H2, label='H2[v] in Fulcher (to be compared against experiment)')
    plt.xlabel('Vibrational state')
    plt.ylabel('H2[v]/H2[0]')
    plt.legend()

    plt.show()



input_crm = 'input_simple.dat'
#simple input file with specified reactions (input_simple is only electron-impact interactions)
run_demo(input_crm)
#run demo with deuterium mass rescaling (e.g. molecular CX depends on relative velocity D+ and D2, which is different (at same Ti) for D than H)
run_demo(input_crm,iso_mass=2)
input_crm = 'input.dat'
#simple input file with additional dissociation & molecular depletion channels
run_demo(input_crm)