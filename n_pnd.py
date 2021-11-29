import numpy as np

def rtc(v): # row vector to column vector
    return v[:,np.newaxis]
def coef_mat(n):
    m = np.zeros((n,n))
    for i in range(n):
        for j in range(n):
            m[i,j] = (n-max(i,j))/(n-i)
    return m
def setup_coef(N):
    c = {}
    for i in range(len(N)):
        c[N[i]] = coef_mat(N[i])
    return c

def theta_mat(t):
    m = np.reshape(np.repeat(t, t.size), (t.size,t.size))
    return m-m.T
def motion_mat(t, coef):
    mat = coef[t.size]*np.cos(theta_mat(t))
    return np.linalg.inv(mat)
def motion_vect(t, o, G, coef):
    v = -G*np.sin(t)
    o_2 = np.square(o)
    th_mat = np.sin(theta_mat(t))
    for i in range(t.size):
        for j in range(t.size):
            delta = coef[t.size][i][j]*th_mat[j][i]*o_2[j]
            v[i] += delta
    return v
def lagrange_f(t, o, G, coef):
    return [o, np.matmul(motion_mat(t, coef), rtc(motion_vect(t, o, G, coef))).flatten()]


def evolve_pnd(t0, o0, T, dt, G, coef):
    theta = t0
    omega = o0
    data_x = []
    data_y = []
    data_t = []
    data_o = []
    for i in range(int(T/dt)):
        data_x.append(np.append(0, np.cumsum(np.sin(theta))))
        data_y.append(np.append(0, np.cumsum(-np.cos(theta))))
        data_t.append(np.copy(theta))
        data_o.append(np.copy(omega))
        k1 = lagrange_f(theta, omega, G, coef)
        k2 = lagrange_f(theta + 0.5*dt*k1[0], omega+0.5*dt*k1[1], G, coef)
        k3 = lagrange_f(theta + 0.5*dt*k2[0], omega+0.5*dt*k2[1], G, coef)
        k4 = lagrange_f(theta + dt*k3[0], omega+dt*k3[1], G, coef)
        d_theta = k1[0]+2*k2[0]+2*k3[0]+k4[0]
        d_omega = k1[1]+2*k2[1]+2*k3[1]+k4[1]
        theta += d_theta*dt/6.0
        omega += d_omega*dt/6.0
    return [data_x, data_y, data_t, data_o]

def sep_vect(data_t_0, data_t_1): # separation vector between two pendulums over time
    data_t_0 = np.asarray(data_t_0)
    data_t_1 = np.asarray(data_t_1)
    data_t_0 = np.append(data_t_0[0], data_t_0[1], axis=1)
    data_t_1 = np.append(data_t_1[0], data_t_1[1], axis=1)
    z = np.linalg.norm(np.abs(data_t_1 - data_t_0), axis=1)
    return z
def lyapunov(data_t_0, data_t_1, dt): # calculations associated with the Lyapunov exponent
    z = sep_vect(data_t_0, data_t_1)
    log_z = np.log(z) - np.log(z[0])
    time = np.arange(log_z.size)*dt
    l_approx = np.divide(log_z[1:], time[1:])
    # take limit for long time and small initial separation
    l_exp = np.mean(l_approx[int(-1.0/dt):]) # take mean of Lyapunov exponent calculations in the last 1 s to reduce noise
    return [time[1:], l_approx, l_exp]
