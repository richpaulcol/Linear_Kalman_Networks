import numpy as np
import filterpy as fp
import pylab as pp
from filterpy.kalman import KalmanFilter
import pykalman as pk

import scipy.signal as ss

#####	Very Simple Model to explore
X0 = np.array([1,1])
A = np.array([[0.99,0.5],[-0.1,0.95]])

P = np.array([[1,0],[0,1]])
Q = np.array([[0.5,0.1],[0.1,0.7]])

Q = Q


#####  	Generation of some data
#Noise vector
QL = np.linalg.cholesky(Q)
PL = np.linalg.cholesky(P)

End = 100
np.random.seed(3)

X = np.zeros((End,2))
X[0] = np.dot(PL,np.random.normal(0,1,2)) + X0

for i in range(1,End):
	w = np.random.normal(0,1,2)
	w = np.dot(QL,w)
	X[i] = np.dot(A,X[i-1]) + w
	
X_measured = X[:,0] + np.random.normal(0,0.2,End)

##### A Kalman Filter using that data
kf = KalmanFilter(dim_x=2, dim_z=1)
kf.x = X0
kf.F = A
kf.H = np.array([[1.,0.]])
kf.R = np.array([[20.]])
kf.alpha = 1.

kf.P = P
kf.Q = Q

X_kf = np.zeros((End,2))
X_kf[0] = kf.x
for i in range(1,End):
	kf.predict()
	kf.update(X_measured[i])
	X_kf[i] = kf.x





##### pyKalman implementation
Ap = A+np.random.normal(0,0.01,[2,2])
kf2 = pk.KalmanFilter(transition_matrices=Ap,observation_matrices=kf.H, transition_covariance=Q,observation_covariance=kf.R, initial_state_mean=X0,initial_state_covariance=P)
X_pkf = np.zeros((End,2))
X_pkf2 = np.zeros((End,2))
X_pkf[0] = X0
X_pkf2[0] = X0

X_pkf[:]= kf2.filter(X_measured)[0]


kf2.em_vars=['transition_matrices']#, 'initial_state_mean', 'initial_state_covariance']

kf2.em(X[:,0],n_iter = 10)
X_pkf2[:]= kf2.filter(X_measured)[0]

f,axs = pp.subplots(2,1,sharex=True)
for i in range(2):
	axs[i].plot(X[:,i],label = 'Data')
	axs[i].plot(X_kf[:,i],label = 'Filter')
	axs[i].plot(X_pkf[:,i],label = 'pyKalman Filter')
	axs[i].plot(X_pkf2[:,i],label = 'EM pyKalman Filter')
axs[0].plot(X_measured,label ='Measured Data')
axs[0].legend()

pp.show()

