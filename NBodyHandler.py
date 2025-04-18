
# import relevant packages
import numpy as np
from scipy.integrate import solve_ivp


# define class
class NBodyHandler:
    """
    This class is used to handle the n-bodies, store their positions,
    calculate each of their accelerations on each other, and finally
    integrate the motion of each of the bodies.
    """

    def __init__(self, mass_array, position_array, velocity_array,
                 G=1, epsilon=0):
        """
        Inputs:
            mass_array: np.array of size N (the number of bodies)
            position_array: np.array of size (N,2) with initial positions
            velocity_array: np.array of size (N,2) with initial velocities

        Optional inputs:
            G: gravitational constant. default
            epsilon: small offset between all bodies which prevents
                singularities. default=0
        """
        self.mass = np.array(mass_array)
        self.position = np.array(position_array)
        self.velocity = np.array(velocity_array)
        self.G = G
        self.N = len(self.mass)
        self.epsilon = epsilon

        # create the relevant flattened arrays
        self.fposition = np.zeros(2*self.N)
        self.fmomentum = np.zeros(2*self.N)
        for idx in range(0, self.N):
            self.fposition[2*idx] = self.position[idx][0]
            self.fposition[2*idx+1] = self.position[idx][1]

            self.fmomentum[2*idx] = self.velocity[idx][0]*self.mass[idx]
            self.fmomentum[2*idx+1] = self.velocity[idx][1]*self.mass[idx]

        # this is the concat of position and momentum - used when integrating
        self.fvec = np.concat((self.fposition, self.fmomentum))

    def get_rhs(self, t, vec):
        """
        Computes the right-hand side of the equations of motion in the
        Hamiltonian formalism (see ipynb for derivation).

        Inputs:
            t: time (but we have no time dependence)
            vec: the np.array of size 4N containing values like:
                vec = [x1, y1, ... xN, yN, px1, py1, ... pN]

        Outputs:
            rhs: np.array of size 4N containing the time derivatives of vec
        """
        # We proceed by computing the first and second halves of the rhs, then
        # taking a concat of them. We first compute the first half of the rhs.
        rhs1 = np.zeros(2*self.N)
        for idx in range(0, self.N):
            rhs1[2*idx] = vec[2*self.N + 2*idx]/self.mass[idx]
            rhs1[2*idx+1] = vec[2*self.N + 2*idx+1]/self.mass[idx]

        # then we compute the second half of the rhs
        rhs2 = np.zeros(2*self.N)
        for idx in range(0, self.N):

            # we need to sum over all of the interactions
            # sumx are the x components, vice versa for sumy
            sumx = 0
            sumy = 0
            for idy in range(0, self.N):
                if idx != idy:
                    sumx += self.mass[idy]*(vec[2*idx] - vec[2*idy])    \
                        * np.pow((vec[2*idx] - vec[2*idy])**2
                                 + (vec[2*idx + 1] - vec[2*idy + 1])**2
                                 + self.epsilon**2, -3/2)
                    sumy += self.mass[idy]*(vec[2*idx+1] - vec[2*idy+1])\
                        * np.pow((vec[2*idx] - vec[2*idy])**2
                                 + (vec[2*idx + 1] - vec[2*idy + 1])**2
                                 + self.epsilon**2, -3/2)

            # here we set the x and y variables
            rhs2[2*idx] = -self.G*self.mass[idx]*sumx
            rhs2[2*idx+1] = -self.G*self.mass[idx]*sumy

        # now we concat rhs 1 and 2 and return them
        return np.concat((rhs1, rhs2))

    def hamiltonian(self, fvec):
        """
        Computes the Hamiltonian for a given flattened fvec.

        Inputs:
            fvec: a flattened vec -- see __init__ for what this means

        Outputs:
            H: the hamiltonian
        """
        # sum up kinetic energies
        T = 0
        for idx in range(0, self.N):
            T += (fvec[2*self.N + 2*idx]**2 + fvec[2*self.N + 2*idx+1]**2) \
                / (2*self.mass[idx])

        # sum up potential energies
        U = 0
        for idx in range(0, self.N):
            for idy in range(idx+1, self.N):
                U += -self.G*self.mass[idx]*self.mass[idy] \
                    / np.sqrt((fvec[2*idx] - fvec[2*idy])**2
                              + (fvec[2*idx+1] - fvec[2*idy+1])**2
                              + self.epsilon**2)

        H = T + U
        return H

    def euler(self, func, t_arr):
        """
        Implemented the Euler method for integrating motion.

        Inputs:
            func: A function returning the rhs of the differential equation
                d/dt y = func
            t_arr: An array of time-values to integrate over. It is best if
                all values are equally spaced, but it is not necessary.

        Output:
            t_arr: the original t array
            vec: An array of size(4*N, len(t_arr)) containing the integrated
                solution
        """
        # initialize vec array
        vec = np.zeros((4*self.N, len(t_arr)))
        vec[:, 0] = self.fvec

        for idx in range(1, len(t_arr)):
            # calculate time step
            dt = t_arr[idx] - t_arr[idx-1]

            # find vec at next time step
            vec[:, idx] = vec[:, idx-1] + dt*func(t_arr[idx-1], vec[:, idx-1])

        return t_arr, vec

    def leapfrog(self, func, t_arr):
        """
        Implemented the leapfrog method for integrating motion.
        Note: this implementation is really inefficient because it computes
        the potential energy more often than it needs to.
        https://en.wikipedia.org/wiki/Leapfrog_integration

        Inputs:
            func: A function returning the rhs of the differential equation
                d/dt y = func
            t_arr: An array of time-values to integrate over. It is best if
                all values are equally spaced, but it is not necessary.

        Output:
            t_arr: the original t array
            vec: An array of size(4*N, len(t_arr)) containing the integrated
                solution
        """
        # initialize vec array
        vec = np.zeros((4*self.N, len(t_arr)))
        vec[:, 0] = self.fvec

        for idx in range(1, len(t_arr)):
            # calculate time step
            dt = t_arr[idx] - t_arr[idx-1]

            # get pos and mom arrays
            pos_arr = vec[:2*self.N, idx-1]
            mom_arr = vec[2*self.N:, idx-1]

            # get half step of momentum
            mom_half_arr = mom_arr + 0.5*func(t_arr[idx-1],
                                              vec[:, idx-1])[2*self.N:]*dt

            # find vec at next time step -- update pos, mom arrays
            vec_mom_arr = np.concat((pos_arr, mom_half_arr))
            pos_arr = pos_arr + func(t_arr[idx-1],
                                     vec_mom_arr)[:2*self.N]*dt

            vec_pos_arr = np.concat((pos_arr, mom_half_arr))
            mom_arr = mom_half_arr + 0.5*func(t_arr[idx-1],
                                              vec_pos_arr)[2*self.N:]*dt

            # save to vec as concat of both pos, mom arrs
            vec[:, idx] = np.concat((pos_arr, mom_arr))

        return t_arr, vec

    def leapfrog4(self, func, t_arr):
        """
        Implemented the 4th order leapfrog method using Yoshida coefficients.
        Note: this implementation is really inefficient because it computes
        the potential energy more often than it needs to.
        https://en.wikipedia.org/wiki/Leapfrog_integration#4th_order_Yoshida_integrator

        Inputs:
            func: A function returning the rhs of the differential equation
                d/dt y = func
            t_arr: An array of time-values to integrate over. It is best if
                all values are equally spaced, but it is not necessary.

        Output:
            t_arr: the original t array
            vec: An array of size(4*N, len(t_arr)) containing the integrated
                solution
        """
        # initialize vec array
        vec = np.zeros((4*self.N, len(t_arr)))
        vec[:, 0] = self.fvec

        # calculate yoshida coefficients
        w0 = -np.pow(2, 1/3)/(2 - np.pow(2, 1/3))
        w1 = 1/(2 - np.pow(2, 1/3))
        c1 = w1/2
        c4 = c1
        c2 = (w0 + w1)/2
        c3 = c2
        d1 = w1
        d3 = d1
        d2 = w0

        for idx in range(1, len(t_arr)):
            # calculate time step
            dt = t_arr[idx] - t_arr[idx-1]

            # get pos and mom arrays
            x0 = vec[:2*self.N, idx-1]
            p0 = vec[2*self.N:, idx-1]

            # first step: x1 and p1
            x1 = x0 + c1*func(t_arr[idx-1],
                              vec[:, idx-1])[:2*self.N]*dt
            vec1 = np.concat((x1, p0))
            p1 = p0 + d1*func(t_arr[idx-1],
                              vec1)[2*self.N:]*dt
            vec1 = np.concat((x1, p1))

            # second step: x2 and p2
            x2 = x1 + c2*func(t_arr[idx-1],
                              vec1)[:2*self.N]*dt
            vec2 = np.concat((x2, p1))
            p2 = p1 + d2*func(t_arr[idx-1],
                              vec2)[2*self.N:]*dt
            vec2 = np.concat((x2, p2))

            # third step: x3, p3
            x3 = x2 + c3*func(t_arr[idx-1],
                              vec2)[:2*self.N]*dt
            vec3 = np.concat((x3, p2))
            p3 = p2 + d3*func(t_arr[idx-1],
                              vec3)[2*self.N:]*dt
            vec3 = np.concat((x3, p3))

            # final step
            x4 = x3 + c4*func(t_arr[idx-1],
                              vec3)[:2*self.N]*dt
            p4 = p3
            vec[:, idx] = np.concat((x4, p4))

        return t_arr, vec

    def solve(self, t_eval,
              method="RK45", rtol=1e-3, atol=1e-6):
        """
        Integrates the equation of motion over some time span.

        Inputs:
            t_span: the time values over which to integrate
            method: the kind of integrator to use. options are
                RK45: use the standard solve_ivp method
                Euler: use this program's implementation of Euler
                Leapfrog: use this program's implementation of leapfrog
                Leapfrog4: use this program's implementation of fourth
                    order leapfrog using yoshida coefficients

        Outputs:
            t: The values of t used to integrate
            vec: the noninterpolated values of vec

        Additionally, this variable is set:
            self.sol: an OdeSolution object which returns an interpolated
                solution -- only created when method=RK45
        """

        # pick an integrator
        if method == "RK45":
            returned_solution = solve_ivp(self.get_rhs,
                                          (t_eval[0], t_eval[-1]),
                                          self.fvec, method="RK45",
                                          t_eval=t_eval,
                                          dense_output=True,
                                          rtol=rtol, atol=atol)
            self.t = returned_solution.t
            self.vec = returned_solution.y
            self.sol = returned_solution.sol

        elif method == "Euler":
            self.t, self.vec = self.euler(self.get_rhs, t_eval)

        elif method == "Leapfrog":
            self.t, self.vec = self.leapfrog(self.get_rhs, t_eval)

        elif method == "Leapfrog4":
            self.t, self.vec = self.leapfrog4(self.get_rhs, t_eval)

        else:
            print("Error: invalid integration method. defaulting to leapfrog.")
            self.t, self.vec = self.leapfrog(self.get_rhs, t_eval)

        # compute hamiltonian for the solution
        self.ham = np.zeros(len(t_eval))
        for idx in range(0, len(t_eval)):
            self.ham[idx] = self.hamiltonian(self.vec[:, idx])

        return self.t, self.vec
