__author__ = 'mpopovic'

from lib import *

import numpy as np
import matplotlib.pyplot as plt
import scipy.integrate as it

## READING EXPERIMENTAL DATA
e_time, e_shear, e_elon = [], [], []
inp = open('WT.txt', 'r')
for line in inp.readlines():
    dat = line.rstrip().split(',')
    e_time.append(float(dat[0]))
    e_shear.append(float(dat[1]))
    e_elon.append(float(dat[2]))
inp.close()
e_time = [x-16.0 for x in e_time]

e_area, e_number = [], []
inp = open('area_and_cell_number.dat','r')
for line in inp.readlines():
    dat = line.rstrip().split(',')
    e_area.append(float(dat[0]))
    e_number.append(float(dat[1]))
inp.close()

e_relative_area = np.array(e_area[6:-5])/e_area[0]
e_relative_number = np.array(e_number[6:-5])/e_number[0]
e_relative_area_times_number = e_relative_area*e_relative_number

##PARAMETERS
L, Lx = 0.38, 1.
Gamma =  200
tau, tauT1 = 1.9, 3.9
tauA_bar = 0.8
t0_bar = 1.8
k1, k2 = 2 , 4
zetaBar = 5
c = -0.17
xi = c/tau
xi_hinge = -0.22/tau
tauA_xi = 0.05
t0_xi = -0.3
xi_undelayed_frac = 0.1

simple_wing_parameters = Wing_parameters(Gamma, tau, tauT1,k1,k2,zetaBar,xi, xi_undelayed_frac, xi_hinge,tauA_bar,t0_bar,tauA_xi,t0_xi,L,Lx)

## VXX LOGISTIC FUNCTION
# time_range_vxx = np.arange(0.1, 19, 0.1)
# shear = 0.5*integrate_v_blade_basic(time_range_vxx, simple_wing_parameters)
# compare_plot([time_range_vxx, e_time], [shear, e_shear], simple_wing_parameters, name = 'shear', x_label='time[h] - 16', y_label='shear', save_path='./shear/', save = True)

time_range_vxx = np.arange(0.1, 19, 0.1)
shear = 0.5*integrate_v_blade(time_range_vxx, simple_wing_parameters)
compare_plot([time_range_vxx, e_time], [shear, e_shear], simple_wing_parameters, name = 'shear', x_label='time[h] - 16', y_label='shear', save_path='./shear/')



## Qm LOGISTIC FUNCTION
# time_range = np.arange(0.1, 19, 0.1)
# deformation = 0.5*integrate_Qm_blade_basic(time_range, simple_wing_parameters)
# compare_plot([time_range, e_time], [deformation, e_elon], simple_wing_parameters, name = 'elongation', x_label='time[h] - 16', y_label='elongation', save_path='./deformation/', save = True)

time_range = np.arange(0.1, 19, 0.1)
deformation = 0.5*integrate_Qm_blade(time_range, simple_wing_parameters)
compare_plot([time_range, e_time], [deformation, e_elon], simple_wing_parameters, name = 'elongation', x_label='time[h] - 16', y_label='elongation', save_path='./deformation/')

## Qp LOGISTIC FUNCTION
time_range_qp = e_time
Qp = integrate_Qp_blade_basic(time_range_qp, simple_wing_parameters)
compare_plot([e_time, e_time, e_time], [np.exp(Qp), np.exp(Qp)/e_relative_number, e_relative_area] ,simple_wing_parameters, label = ['1D theory', 'theory with CD correction', 'experiment'], name = 'area', x_label='time[h] - 16', y_label='A/A0', save_path='./area/', save = True)


## ONLY STRESS DELAY
deformation_list = []
time_range_list = []
label_list = []
time_range = np.arange(0.1,19,0.1)
for xi_undelayed_frac in np.arange(0.0, 1.2, 0.2):
    print xi_undelayed_frac
    simple_wing_parameters = Wing_parameters(Gamma, tau, tauT1,k1,k2,zetaBar,xi, xi_undelayed_frac, xi_hinge,tauA_bar,t0_bar,tauA_xi,t0_xi,L,Lx)
    time_range_list.append(time_range)
    deformation_list.append(0.5*integrate_Qm_blade(time_range, simple_wing_parameters))
    label_list.append(str(xi_undelayed_frac))
deformation_list.append(e_elon)
time_range_list.append(e_time)
label_list.append('experiment')
compare_plot(time_range_list, deformation_list, simple_wing_parameters, label = label_list)

## ONLY T1 DELAY
deformation_list = []
time_range_list = []
label_list = []
time_range = np.arange(0.1,19,0.1)
for tauT1 in np.arange(3.6, 4.6, 0.2):
    print tauT1
    simple_wing_parameters = Wing_parameters(Gamma, tau, tauT1,k1,k2,zetaBar,xi, xi_undelayed_frac, xi_hinge,tauA_bar,t0_bar,tauA_xi,t0_xi,L,Lx)
    time_range_list.append(time_range)
    deformation_list.append(0.5*integrate_Qm_blade(time_range, simple_wing_parameters))
    label_list.append(str(tauT1))
deformation_list.append(e_elon)
time_range_list.append(e_time)
label_list.append('experiment')
compare_plot(time_range_list, deformation_list, simple_wing_parameters, label = label_list)


## ONLY T1 DELAY AND STRESS DELAY

time_range = np.arange(0.1,19,0.1)
for xi_undelayed_frac in np.arange(0.0, 1.2,0.2):
    deformation_list = []
    time_range_list = []
    label_list = []
    for tauT1 in np.arange(3.6, 4.6, 0.2):
        print xi_undelayed_frac, tauT1
        simple_wing_parameters = Wing_parameters(Gamma, tau, tauT1,k1,k2,zetaBar,xi, xi_undelayed_frac, xi_hinge,tauA_bar,t0_bar,tauA_xi,t0_xi,L,Lx)
        time_range_list.append(time_range)
        deformation_list.append(0.5*integrate_Qm_blade(time_range, simple_wing_parameters))
        label_list.append(str(tauT1))
    deformation_list.append(e_elon)
    time_range_list.append(e_time)
    label_list.append('experiment')
    compare_plot(time_range_list, deformation_list, simple_wing_parameters, label = label_list, name='undelayed_frac_'+str(xi_undelayed_frac), save_path='./', save = True, show = False)



## CALCULATING VXX in hinge
time_range_vxx = np.arange(0.1, 19, 0.1)
result_vxx_hinge = [it.quad(lambda x: 1./2./np.pi*np.real(calculate_v_hinge(x, Gamma, tau, tauT1, tauA,  k1, k2, zetaBar) * np.exp(np.complex128(1.j)*x *(t_vxx-t0))), 0, 100) for t_vxx in time_range_vxx]
plt.figure()
plt.plot(time_range_vxx, zip(*result_vxx_hinge)[0], label = '1D theory')
plt.plot(e_time, e_shear, label= 'experiment')
plt.legend()
plt.xlabel('time[h] - 16', fontsize=20)
plt.ylabel('half_v_xx', fontsize=20)
plt.show()
# plt.savefig('half_vxx_logistic_L_'+str(L)+'_Gamma_'+ str(Gamma)+'_tau_'+str(tau)+'_tauT1_'+str(tauT1)+'_k1_'+str(k1)+'_k2_'+str(k2)+'_zetaBar_'+str(zetaBar)+'.png')

w_range = np.arange(0.0001, 1, 0.0001)
plt.figure()
plt.plot(w_range, calculate_v_blade(w_range, simple_wing_parameters))
plt.xscale('log')
plt.show()

w_range = np.arange(0.01, 10, 0.01)
plt.figure()
plt.plot(w_range, np.abs(calculate_Qm_vPart_blade(w_range, simple_wing_parameters)))
plt.yscale('log')
plt.show()

##CHECK SELF CONSISTENCY
Qm = np.array(zip(*result)[0])
Q_deriv = 10.*(Qm[1:]-Qm[:-1])
plt.figure()
plt.plot(time_range, Qm/tau)
plt.plot(time_range, zip(*last_term)[0])
plt.show()


plt.figure()
plt.plot(time_range, zip(*(result))[0], label = 'total')
plt.plot(e_time, e_elon, label = 'experiment')
plt.plot(time_range, zip(*(result-last_term))[0], label = 'isotropic active contribution')
#plt.plot(time_range, zip(*last_term)[0])
plt.plot(time_range, zip(*(last_term))[0], label = 'anisotropic active contribution')

# plt.plot(time_range, zip(*(result_noT1_delay-last_term_nodelay))[0], label = 'without T1 delay')
# plt.plot(time_range, zip(*last_term_nodelay)[0])
plt.legend(fontsize = 11)
plt.ylim(-0.05, 0.35)
plt.xlabel('time', fontsize=20)
plt.ylabel('Q_m', fontsize=20)
plt.savefig('half_Qm_logistic_contributions_withoutT1_L_'+str(L)+'_Gamma_'+ str(Gamma)+'_tau_'+str(tau)+'_tauT1_'+str(tauT1)+'_k1_'+str(k1)+'_k2_'+str(k2)+'_zetaBar_'+str(zetaBar)+'.png')
plt.show()

av = np.mean[it.quad(lambda x: np.real(calculate_v_step(x, Gamma, tau, tauT1,  k1, k2) * np.exp(np.complex128(1.j)*x *(t))), 0, 1000+y) for y in np.arange(0,100,1)]
np.mean(av)
plt.figure()
plt.plot(range(len(av)), av)
plt.show()



## STEP FUNCTION
time_range = np.arange(0.1, 15, 0.1)
result = [it.quad(lambda x: np.real(calculate_v_logistic(x, Gamma, tau, tauT1,  k1, k2) * np.exp(np.complex128(1.j)*x *(t))), 0, 100) for t in time_range]
result_noT1_delay = [it.quad(lambda x: np.real(calculate_v_logistic(x, Gamma, tau, 0.,  k1, k2) * np.exp(np.complex128(1.j)*x *(t))), 0, 100) for t in time_range]
n_smooth = 1
kernel = np.ones(n_smooth)/n_smooth
plt.figure()
plt.plot(time_range, np.convolve(zip(*result)[0], kernel,'same'), label = 'with T1 delay')
plt.plot(time_range, np.convolve(zip(*result_noT1_delay)[0], kernel,'same'), label = 'without T1 delay')
plt.legend()
plt.xlabel('time', fontsize=20)
plt.ylabel('v_xx', fontsize=20)
plt.savefig('vxx_stepish_logistic.png')
plt.show()


## TESTING Qm -logistic
w_range = np.arange(0.01,10,0.01)
w_minus_range = np.arange(-1, 0., 0.01)
Qm_logistic_test = calculate_Qm_logistic(w_range, Gamma, tau, tauT1, tauA, k1, k2, zetaBar)
Qm_logistic_test_minus = calculate_Qm_logistic(w_minus_range, Gamma, tau, tauT1, tauA, k1, k2, zetaBar)
Qm_logistic_test-Qm_logistic_test_minus
plt.figure()
plt.plot(w_range, np.real(Qm_logistic_test))
plt.plot(w_minus_range, np.real(Qm_logistic_test_minus))
plt.show()

## TESTING V(\omega) - logistic activation
t=0
np.real(it.quad(lambda x: (calculate_v_logistic(x, Gamma, tau, tauT1, tauA,  k1, k2, zetaBar)-calculate_v_logistic(-x, Gamma, tau, tauT1, tauA,  k1, k2, zetaBar)) * np.exp(np.complex128(1.j)*x *(t-t0)), 0, 1))
it.quad(lambda x: np.real((calculate_v_logistic(x, Gamma, tau, tauT1, tauA,  k1, k2, zetaBar)-calculate_v_logistic(-x, Gamma, tau, tauT1, tauA,  k1, k2, zetaBar)) * np.exp(np.complex128(1.j)*x *(t-t0))), -1, 0)
w_range = np.arange(-1,1,0.01)
v_logistic_test = calculate_v_logistic(w_range, Gamma, tau, tauT1, tauA, k1, k2, zetaBar)* np.exp(np.complex128(1.j)*w_range *(t-t0))
v_logistic_test_inverse = calculate_v_logistic(-w_range, Gamma, tau, tauT1, tauA, k1, k2, zetaBar)* np.exp(np.complex128(1.j)*w_range *(t-t0))
plt.figure()
plt.plot(w_range, np.real(v_logistic_test))
plt.plot(w_range, np.real(v_logistic_test_inverse))
#plt.plot(w_range, v_logistic_test[::-1])
#plt.plot(w_range, np.imag(v_logistic_test))
plt.show()

## TESTING V(\omega) - step activation
w_range = np.logspace(-2,5,num=100, endpoint= True)
v_step_test = calculate_v_step(w_range, Gamma, tau, tauT1, k1, k2)
plt.figure()
plt.plot(w_range, np.real(v_step_test))
plt.plot(w_range, np.abs(np.imag(v_step_test)))
plt.xscale('log')
plt.yscale('log')
plt.show()


## TESTING LAMBDA
w_range = np.arange(0.01,10,0.01)
lambda_test = Lambda(w_range, Gamma, tau, tauT1, k1, k2)
plt.figure()
plt.plot(w_range, np.real(lambda_test))
plt.plot(w_range, np.imag(lambda_test))
plt.show()



##checking dual margin data
# check_time, check_shear = [], []
# inp = open('data/111102/dualMarginDeformation.dat','r')
# inp.readline()
# for line in inp.readlines():
#     dat = line.rstrip().split()
#     check_time.append(float(dat[0]))
#     check_shear.append(0.5*(float(dat[2])-float(dat[5])))
# inp.close()
# plt.figure()
# plt.plot(check_time, np.convolve(check_shear,np.ones(5)/5, 'same'))
# plt.show()


# coth = lambda x: 1/np.tanh(x)
