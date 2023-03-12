#script for calculating effective rates based on CRUMPET

import CRUMPET
import numpy as np

input_crm = 'input_3.dat' #simple input file with only e- + H2(v1) -> e- + H2(v2), where abs(v1-v2) = 1
indx_H2v = np.append(2,np.arange(3,17))

#initialise CRM
crm = CRUMPET.Crumpet(input_crm)

#make Te & ne vectors
Tev = np.linspace(0.2,10,100)
Tiv = Tev
ne = 1e19*np.ones(np.shape(Tev))
crm.source[2] = 1e-100

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
    vibr_resolved_CX[i,:] = eval_1D(X.reactions['H2VIBR']['H.2']['2.'+str(i)+'L2'],Tiv)

#Now use fv_H2 as a weight and sum the total reaction rate to generate the effective rate
eff_mol_cx = np.sum(vibr_resolved_CX*fv_H2,axis=0)

#compare effective rate against tabulated one
h3_2_3_tab = eval_1D(X.reactions['AMJUEL']['H.2']['3.2.3'],Tev)

import matplotlib.pyplot as plt
plt.figure()
plt.loglog(Tev,np.transpose(fv_H2))

plt.figure()
plt.loglog(Tev,eff_mol_cx,label='Effective rate (NEW)')
plt.loglog(Tev,h3_2_3_tab,label='Tabulated rate')
plt.legend()
plt.show()

print('...')


