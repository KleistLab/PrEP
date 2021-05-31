#!/usr/bin/env python3
# Lanxin Zhang

import numpy as np
from scipy.integrate import solve_ivp


class Sample:
    """class for samples with different PK parameters

    Args:
    k_a: absorption rate constant [1/h]
    V_p: Vp/Fbio, bioavailability-adjusted volume of peripheral compartment [L]
    V_c: Vc/Fbio, bioavailability-adjusted volume of central compartment [L]
    CL: CL/Fbio, bioavailability-adjusted drug clearance [L/h]
    Q: Q/Fbio, bioavailability-adjusted intercompartmental clearance rate [L/h]
    """

    M_DTG = 419.38      # mol mass of DTG

    def __init__(self, k_a, V_p, Q, CL, V_c):
        """initialize the sample object.

        Attributes:
        solutions: store all OdeSolution for the concentration.
        day_array: store the days on those drug is taken (denote as 1, else 0)
        mdose: mass of each dose.
        countdose: count of doses in the regimen.
        Dmax_array: store the maximal concentration and
            the corresponding time point of each day  [(t, Dmax)]
            if drug is taken on this day,  else 0.  (Length = countdose)
        """

        self.__k_a = k_a
        self.__V_p = V_p
        self.__Q = Q
        self.__CL = CL
        self.__V_c = V_c
        self.day_array = list()
        self.solutions = list()
        self.mdose = 0
        self.adh = 0
        self.countdose = 0
        self.Dmax_array = list()

    def reset(self):
        """reset the object variables"""

        self.day_array = list()
        self.solutions = list()
        self.mdose = 0
        self.adh = 0
        self.countdose = 0
        self.Dmax_array = list()

    def __pkModel(self, t, y):

        """PK model ODE used for solving ODE with scipy,
        return a vector of the right hand side of  the ODEs.

        t: time
        y: vector [Z1, Z2, Z3]
        """
        z1, z2, z3 = y[0], y[1], y[2]
        dydt = [-self.__k_a * z1,
                self.__k_a / self.__V_c * z1 - self.__CL / self.__V_c * z2
                - self.__Q / self.__V_c * z2 + self.__Q / self.__V_p * z3,
                self.__Q / self.__V_c * z2 - self.__Q / self.__V_p * z3]
        return dydt

    def __solveODE(self, t_start, t_end, n_dosing, c_plasma, c_peri):

        """solve ODE for drug concentration based on the PK parameters of
        a sample, return bunch object containing solutions.

        t_start, t_end: time span for the integration
        n_dosing: amount of oral taken DTG [nmol]
        c_plasma: concentration of DTG in blood [nM]
        c_peri: concentration of DTG in peripheral compartment [nM]

        """
        t = [t_start, t_end]
        y0 = [n_dosing, c_plasma, c_peri]
        k_a = self.__k_a
        V_p = self.__V_p
        Q = self.__Q
        CL = self.__CL
        V_c = self.__V_c
        sol = solve_ivp(self.__pkModel, t, y0, dense_output=True)
        return sol

    def getConcentration(self, m_DTG, tp, adherence, count_doses, adh_pattern=None, day_upperbound=100):

        """pharmacokinetic profile of DTG, store the ODESolution objects of
         each day in self.solutions. If no drug is taken for today, the
         ODE solution object will be inherited form the last day, but the
         time span for integration will be extended; if no drug is taken
         for all days before, None will be stored in self.solutions.

         return value: concentration of DTG at time point tp

         m_DTG: mass of oral taken DTG [ng]
         tp: time point relative to first dosing [hr] (exposure occurs)
            (can be negative, time before dosing)
         adherence: possibility to take drug
         count_doses: count of doses in a regimen
         adh_pattern: a list contains 0 and 1 for adherence pattern, length is count_doses. 
             if adh_pattern is given, adherence is ignored.
         day_upperbound: upperbound of days, calculate until the count_doses+day_upperbound 
             is reached

         """

        # day_max: maximum days for simulation, must be larger than count_doses
        day_max = count_doses + day_upperbound

        if tp < 0:
            return 0

        # recalculate the concentrations if one of the inputs is changed
        if not self.mdose or self.mdose != m_DTG or self.adh != adherence\
           or self.countdose != count_doses:
            self.reset()
            self.mdose = m_DTG
            self.adh = adherence
            self.countdose = count_doses
            if adh_pattern is not None:
                self.day_array = adh_pattern
            else:
                # precalculate the days on which the drug is taken according
                # to the adherence, drug must be taken on the first day
                self.day_array = np.random.choice(2, size=count_doses, p=[1-adherence,adherence])
                self.day_array[0] = 1 
            self.day_array = list(self.day_array) + [0] * (day_max - count_doses)

            # compute the ODE solution for the first day
            n_dosing, c_plasma, c_peri = m_DTG / self.M_DTG, 0, 0
            solution_0 = self.__solveODE(0, 24, n_dosing, c_plasma, c_peri).sol
            self.solutions.append(solution_0)
            day = 1
            while day < day_max:
                if self.day_array[day]:
                    n_dosing, c_plasma, c_peri = self.solutions[-1](day*24)[:3]
                    n_dosing += self.mdose / self.M_DTG
                    solution_new = self.__solveODE(day * 24,
                                                   (day + 1) * 24,
                                                   n_dosing,
                                                   c_plasma,
                                                   c_peri).sol
                    self.solutions.append(solution_new)
                    day += 1
                else:
                    beginday = day
                    nodrugdays = 0
                    while day < day_max and not self.day_array[day]:
                        nodrugdays += 1
                        day += 1
                    n_dosing, c_plasma, c_peri = self.solutions[-1](beginday*24)[:3]
                    solution_new = self.__solveODE(beginday * 24,
                                                   day * 24,
                                                   n_dosing,
                                                   c_plasma,
                                                   c_peri).sol
                    self.solutions += [solution_new] * nodrugdays

        daycount = int(tp // 24)
        if daycount >= day_max:
            return 0
        else:
            return self.solutions[daycount](tp)[1]

    def getMaximum(self, m_DTG, tp, adherence, count_doses, threshold=1e-3):
        """compute the maximal concentration of DTG from time point tp on.
        Input variables should be same as those in self.getConcentration
        inside a run.

        return value:  maximal concentration of DTG from time point tp on

        m_DTG: mass of oral taken DTG [ng]
        tp: time point relative to first dosing [hr] (exposure occurs)
            (can be negative, time before dosing)
        adherence: possibility to take drug
        count_doses: count of doses in a regimen
        threshold: used to find maximum [hr]
        """

        def findMaximum(sol, t1, t2):
            """function to find local maximum in a time span using binary search
            return the time point on which the maximum occurs and the maximum

            sol: OdeSolution object
            t1, t2: initial time span for binary search
            """
            t = 0  # time point for the maximum
            # stop until the difference is smaller than threshold
            while t2 - t1 > threshold:
                l, r, m = t1, t2, (t1 + t2) / 2
                l1, r1 = (l + m) / 2, (m + r) / 2
                p = [l, l1, m, r1, r]
                s = list(map(lambda x: sol(x)[1], p))
                if max(s) == sol(m)[1]:
                    t1, t2, t = l1, r1, m
                elif max(s) == sol(l1)[1]:
                    t2, t = m, l1
                elif max(s) == sol(r1)[1]:
                    t1, t = m, r1
                elif max(s) == sol(l)[1]:
                    t2, t = l1, l
                elif max(s) == sol(r)[1]:
                    t1, t = r1, r
            return t, sol(t)[1]
        
        if tp < 0:
            self.getConcentration(m_DTG, count_doses*24, adherence, count_doses)

        # if getConcentration is not called before, call it
        if not self.day_array:
            self.getConcentration(m_DTG, tp, adherence, count_doses)

        # store local maximums in self.Dmax_array to save repetitive calculation
        if not self.Dmax_array:
            for i in range(count_doses):
                if self.day_array[i]:
                    self.Dmax_array.append(findMaximum(self.solutions[i],
                                                       i * 24, (i + 1) * 24))
                else:
                    self.Dmax_array.append((0, 0))
        self.Dmax_array.append((0, 0))

        daycount = int(tp // 24)
        c_curr = self.solutions[daycount](tp)[1]    # concentration at tp
        if daycount >= count_doses:
            return c_curr
        # maximal concentration from today on
        tmax, cmax = max(self.Dmax_array[daycount:], key=lambda x: x[1])

        # if the maximum exists and tp is before the time point of maximum
        if cmax and tp <= tmax:
            return cmax
        # if maximum doesn't exist, i.e. no doses taken after this day,
        # return the concentration of current time point tp
        # since the concentration is monotonically decreasing.
        elif not cmax:
            return c_curr
        # if maximum exists but tp is after that, search the new maximum
        # from next day on and return the bigger one between new maximum and
        # current concentration.
        else:
            cmax = max(self.Dmax_array[daycount+1:], key=lambda x: x[1])[1]
            return max(cmax, c_curr)
