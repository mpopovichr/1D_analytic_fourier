__author__ = 'mpopovic'

import numpy as np
import scipy.integrate as it
import matplotlib.pyplot as plt

class Wing_parameters():
    def __init__(self, Gamma, tau, tauT1, k1, k2, zetaBar, xi, tauA_bar, t0_bar, tauA_xi, t0_xi, L, Lx):
        self.Gamma = Gamma
        self.tau = tau
        self.tauT1 = tauT1
        self.k1 = k1
        self.k2 = k2
        self.zetaBar = zetaBar
        self.tauA_bar = tauA_bar
        self.t0_bar = t0_bar
        self.tauA_xi = tauA_xi
        self.t0_xi = t0_xi
        self.L = L
        self.Lx = Lx
        self.xi = xi


Imag = lambda x: np.complex128(1.j)*x

def Lambda(w, p):
    return np.sqrt(p.Gamma*((Imag(1)*w + 1/((1+Imag(1)*w*p.tauT1)*p.tau))/(p.k1 + p.k2*(1 - Imag(1)/(w*((1+Imag(1)*w*p.tauT1)*p.tau))))))


## BASIC BLADE:
  ## IN FREQUENCY SPACE:
def calculate_v_blade_basic(w, p):
    return np.exp(-Imag(1)*w*p.t0_bar)/((p.Lx-p.L))*Lambda(w, p)/p.Gamma*p.zetaBar*(-1.*Imag(1))*np.pi*p.tauA_bar/(np.sinh(np.pi*w*p.tauA_bar)) * np.sinh(Lambda(w, p)*p.L) * np.sinh(Lambda(w, p)*(p.Lx-p.L))/np.sinh(Lambda(w, p)*p.Lx)

def calculate_Qm_vPart_blade_basic(w, p):
    return 1./(Imag(1)*w+ 1./((1+Imag(1)*w*p.tauT1)*p.tau))*calculate_v_blade_basic(w, p)

def calculate_Qp_blade_basic(w, Gamma, taut, tauT1, tauA, k1, k2, zetabar):
    return 1./(Imag(1)*w)*calculate_v_logistic(w, Gamma, tau, tauT1, tauA, k1, k2, zetaBar)

  ## INTEGRATING TO TIME:
def integrate_v_blade_basic(time_range, p):
    err = 0.0000001
    upper_limit = 1
    ratio = np.abs(np.real(calculate_v_blade_basic(upper_limit, p))/np.real(calculate_v_blade_basic(0.001, p)))
    while ratio > err:
        upper_limit *= 2
        print upper_limit
        ratio = np.abs(np.real(calculate_v_blade_basic(upper_limit, p))/np.real(calculate_v_blade_basic(0.001, p)))
    result_vxx = [it.quad(lambda x: 1./np.pi*np.real(calculate_v_blade_basic(x, p) * np.exp(Imag(1)*x *t_vxx)), 0, upper_limit) for t_vxx in time_range]
    return np.array(zip(*result_vxx)[0])

def integrate_Qm_blade_basic(time_range, p):
    result_Qm = -p.xi*p.tau - np.real(np.array(zip(*2.*p.xi/np.pi*np.array([it.quad(lambda x: -Imag(1)*np.pi*p.tauA_xi/(np.sinh(np.pi*p.tauA_xi*x))/(Imag(1)*x+1./((1+Imag(1)*x*p.tauT1)*p.tau))/(1+Imag(1)*x*p.tauT1)*np.exp(Imag(1)*x*(time-p.t0_xi)),0,100) for time in time_range]))[0]))
    print 'calculated xi part'
    err = 0.0000001
    upper_limit = 1
    ratio = np.abs(np.real(calculate_Qm_vPart_blade_basic(upper_limit, p))/np.real(calculate_Qm_vPart_blade_basic(0.001, p)))
    while ratio > err:
        upper_limit *= 2
        print upper_limit
        ratio = np.abs(np.real(calculate_Qm_vPart_blade_basic(upper_limit, p))/np.real(calculate_v_blade_basic(0.001, p)))
    result_Qm += np.array(zip(*np.real(np.array([it.quad(lambda x: 1/np.pi*(calculate_Qm_vPart_blade_basic(x, p) * np.exp(Imag(1)*x *t_qm)), 0, upper_limit) for t_qm in time_range])))[0])
    return np.array(result_Qm)

aoeuaeouaeu
############################################

## BASIC HINGE:
def calculate_v_hinge_basic(w, Gamma, tau, tauT1, tauA, k1, k2, zetaBar):
    return 1/(L)*Lambda(w, Gamma, tau, tauT1,  k1, k2)/Gamma*zetaBar*(-1.*Imag(1))*np.pi*tauA/(np.sinh(np.pi*w*tauA)) * np.sinh(Lambda(w, Gamma, tau, tauT1,  k1, k2)*L) * np.sinh(Lambda(w, Gamma, tau, tauT1,  k1, k2)*(Lx-L))/np.sinh(Lambda(w, Gamma, tau, tauT1,  k1, k2)*Lx)


#### PLOTTING

def compare_plot(x_t, val_t, x_e, val_e,p,  x_label= '', y_label= '', t_label = 'theory', e_label = 'experiment', save = False, save_path = '', name = ''):
    plt.figure()
    plt.plot(x_t, val_t, label = t_label)
    plt.plot(x_e, val_e, label= e_label)
    plt.legend()
    plt.xlabel(x_label, fontsize=20)
    plt.ylabel(y_label, fontsize=20)
    if save:
        plt.savefig(save_path+name+'_'+str(p.L)+'_Gamma_'+ str(p.Gamma)+'_tau_'+str(p.tau)+'_tauT1_'+str(p.tauT1)+'_k1_'+str(p.k1)+'_k2_'+str(p.k2)+'_zetaBar_'+str(p.zetaBar)+'.png')
    plt.show()