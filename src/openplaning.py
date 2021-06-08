import numpy as np
from scipy import interpolate, signal
from scipy.special import gamma
import ndmath
import warnings

class PlaningBoat():
    """Prismatic planing craft
    
    ATTRIBUTES
    ----------
    speed (float): m/s, speed
    weight (float): N, weight
    beam (float): m, beam
    lcg (float): m, longitudinal center of gravity, measured from the stern
    vcg (float): m, vertical center of gravity, measured from the keel
    r_g (float): m, radius of gyration
    beta (float): deg, deadrise
    epsilon (float): deg, thrust angle w.r.t. keel (CCW with body-fixed origin at 9 o'clock)
    vT (float): m, thrust vertical distance (from keel, positive up)
    lT (float): m, thrust horizontal distance (from stern, positive forward)
    length (float): m, vessel LOA for seaway behavior estimates (default of None)
    H_sig (float): m, significant wave heigth in an irregular sea state (default of None)
    ahr (float): m, average hull roughness (default of 150*10**-6)
    Lf (float): m, flap chord (default of 0)
    sigma (float): dimensionless, flap span-beam ratio (default of 0)
    delta (float): deg, flap deflection (default of 0)
    l_air (float): m, distance from stern to center of air pressure (default of 0)
    h_air (float): m, height from keel to top of square which bounds the air-drag-inducing shape (default of 0)
    b_air (float): m, transverse width of square which bounds the air-drag-inducing shape (default of 0)
    C_shape (float): dimensionless, area coefficient for air-drag-inducing shape. C_shape = 1 means the air drag reference area is h_air*b_air (default of 0)
    C_D (float): dimensionless, air drag coefficient (default of 0.7)
    rho (float): kg/m^3, water density (default of 1025.87)
    nu (float): m^2/s, water kinematic viscosity (default of 1.19*10**-6)
    rho_air (float): kg/m^3, air density (default of 1.225)
    g (float): m/s^2, gravitational acceleration (default of 9.8066)      
    z_wl (float): m, vertical distance of center of gravity to the calm water line (default of 0)
    tau (float): deg, trim angle (default of 5)
    eta_3 (float): deg, additional heave (default to zero)
    eta_5 (float): deg, additional trim (default to zero)
    wetted_lengths_type (int): 1 (default) = Use Faltinsen 2010 wave rise approximation, 2 = Use Savitsky's '64 approach, 3 = Use Savitsky's '76 approach
    z_max_type (int): 1 = Uses 3rd order polynomial fit (faster, recommended), 2 = uses cubic interpolation (table provided for completeness). This is only used if wetted_lenghts_type == 1.
    L_K (float): m, keel wetted length
    L_C (float): m, chine wetted length
    lambda_W (float): dimensionless, mean wetted-length to beam ratio (L_K+L_C)/(2*beam)
    x_s (float): m, distance from keel/water-line intersection to start of wetted chine
    z_max (float): dimensionless, maximum presssure coordinate coefficient (z_max/Ut)
    hydrodynamic_force (numpy float array): [F_x, F_z, M_cg] N, hydrodynamic force (x, y, rot directions in intertial coordinates)
    skin_friction (numpy float array): [F_x, F_z, M_cg] N, skin friction force
    air_resistance (numpy float array): [F_x, F_z, M_cg] N, air resistance force
    flap_force (numpy float array): [F_x, F_z, M_cg] N, flap resultant force
    thrust_force (numpy float array): [F_x, F_z, M_cg] N, thrust resultant force
    net_force (numpy float array): [F_x, F_z, M_cg] N, net force
    mass_matrix (numpy float 2D array): [[A_33 (kg), A_35 (kg*m/rad)], [A_53 (kg*m), A_55 (kg*m^2/rad)]]
    damping_matrix (numpy float 2D array): [[B_33 (kg/s), B_35 (kg*m/(s*rad))], [B_53 (kg*m/s), B_55 (kg*m**2/(s*rad))]]
    restoring_matrix (numpy float 2D array): [[C_33 (N/m), C_35 (N/rad)], [C_53 (N), C_55 (N*m/rad)]]
    porpoising (list): [[eigenvalue result, est. pitch settling time (s)], [Savitsky chart result, critical trim angle (deg)]]
    avg_impact_acc (numpy float array): [n_cg, n_bow] g's, average impact acceleration at center of gravity and bow
    R_AW (float): N, added resistance in waves
    
    METHODS
    -------
    print_description(sigFigs=7, runAllFunctions=True): Return a formatted description of the vessel.
    get_geo_lengths(runGeoLengths=True): This function outputs the geometric lengths. 
    get_forces(runGeoLengths=True): This function calls all the force functions to update the respective object attributes.
    get_steady_trim(x0=[0, 3], tauLims=[0.5, 35], tolF=10**-6, maxiter=50): This function finds and sets the 'equilibrium' point when the vessel is steadily running in calm water.
    get_eom_matrices(runGeoLengths=True): This function returns the mass, damping, and stiffness matrices following Faltinsen 2010.
    check_porpoising(stepEstimateType=1): This function checks for porpoising.
    get_seaway_behavior(): This function calculates the seaway behavior as stated in Savitsky & Brown '76.
    """
    
    def __init__(self, speed, weight, beam, lcg, vcg, r_g, beta, epsilon, vT, lT, length=None, H_sig=None, ahr=150*10**-6, Lf=0, sigma=0, delta=0, l_air=0, h_air=0, b_air=0, C_shape=0, C_D=0.7, z_wl=0, tau=5, rho=1025.87, nu=1.19*10**-6, rho_air=1.225, g=9.8066, wetted_lengths_type=1, z_max_type=1):
        """Initialize attributes for PlaningBoat
        
        PARAMETERS
        ----------
        speed (float): m/s, speed
        
        weight (float): N, weight
        beam (float): m, beam
        lcg (float): m, longitudinal center of gravity, measured from the stern
        vcg (float): m, vertical center of gravity, measured from the keel
        r_g (float): m, radius of gyration

        beta (float): deg, deadrise
        epsilon (float): deg, thrust angle w.r.t. keel (CCW with body-fixed origin at 9 o'clock)
        vT (float): m, thrust vertical distance (from keel, positive up)
        lT (float): m, thrust horizontal distance (from stern, positive forward)
        
        length (float): m, vessel LOA for seaway behavior estimates (default of None)
        H_sig (float): m, significant wave heigth in an irregular sea state (default of None)
        
        ahr (float): m, average hull roughness (default of 150*10**-6)
        
        Lf (float): m, flap chord (default of 0)
        sigma (float): dimensionless, flap span-beam ratio (default of 0)
        delta (float): deg, flap deflection (default of 0)
        
        l_air (float): m, distance from stern to center of air pressure (default of 0)
        h_air (float): m, height from keel to top of square which bounds the air-drag-inducing shape (default of 0)
        b_air (float): m, transverse width of square which bounds the air-drag-inducing shape (default of 0)
        C_shape (float): dimensionless, area coefficient for air-drag-inducing shape. C_shape = 1 means the air drag reference area is h_air*b_air (default of 0)
        C_D (float): dimensionless, air drag coefficient (default of 0.7)
        
        rho (float): kg/m^3, water density (default of 1025.87)
        nu (float): m^2/s, water kinematic viscosity (default of 1.19*10**-6)
        rho_air (float): kg/m^3, air density (default of 1.225)
        g (float): m/s^2, gravitational acceleration (default of 9.8066)
        
        INITIALIZED
        -----------        
        z_wl (float): m, vertical distance of center of gravity to the calm water line (default of 0)
        tau (float): deg, trim angle (default of 5)
        eta_3 (float): deg, additional heave (default to zero)
        eta_5 (float): deg, additional trim (default to zero)
        
        wetted_lengths_type (int): 1 (default) = Use Faltinsen 2010 wave rise approximation, 2 = Use Savitsky's '64 approach, 3 = Use Savitsky's '76 approach
        z_max_type (int): 1 = Uses 3rd order polynomial fit (faster, recommended), 2 = uses cubic interpolation (table provided for completeness). This is only used if wetted_lenghts_type == 1.
        """
        self.speed = speed
        self.weight = weight
        self.beam = beam
        self.lcg = lcg 
        self.vcg = vcg 
        self.r_g = r_g 
        self.beta = beta
        self.epsilon = epsilon 
        self.vT = vT 
        self.lT = lT
        self.length = length
        self.H_sig = H_sig
        self.ahr = ahr
        self.Lf = Lf
        self.sigma = sigma
        self.delta = delta
        self.l_air = l_air
        self.h_air = h_air
        self.b_air= b_air
        self.C_shape = C_shape
        self.z_wl = z_wl
        self.tau = tau
        self.eta_3 = 0
        self.eta_5 = 0
        self.rho = rho
        self.nu = nu
        self.rho_air = rho_air
        self.C_D = C_D
        self.g = g
        
        self.gravity_force = np.array([0, -self.weight, 0])
        
        self.LC_type = wetted_lengths_type
        self.zmax_type = z_max_type
        
    def print_description(self, sigFigs=7, runAllFunctions=True):
        """Return a formatted description of the vessel.
        
        PARAMETERS
        ----------
        sigFigs (integer): Number of significant figures to display (default of 7)
        runAllFunctions (boolean): Runs all functions with default values before printing results.
        """
        if runAllFunctions:
            self.get_geo_lengths()
            self.get_forces(runGeoLengths=False)
            self.get_eom_matrices(runGeoLengths=False)
            self.get_seaway_behavior()
            self.check_porpoising()
        
        volume = self.weight/(self.g*self.rho)
        
        table = [
            ['---VESSEL---'],
            ['Speed', self.speed, 'm/s'],
            ['V_k', self.speed*1.944, 'knot'],
            ['Fn (beam)', self.speed/np.sqrt(self.g*self.beam), ''],
            ['Fn (volume)', self.speed/np.sqrt(self.g*(self.weight/(self.g*self.rho))**(1/3)), ''],
            [''],
            ['Weight', self.weight, 'N'],
            ['Mass', self.weight/self.g, 'kg'],
            ['Volume', self.weight/(self.g*self.rho), 'm\u00B3'],
            ['Beam', self.beam, 'm'],
            ['LCG', self.lcg, 'm from stern'],
            ['VCG', self.vcg, 'm from keel'],
            ['R_g', self.r_g, 'm'],
            ['Deadrise', self.beta, 'deg'], #'\N{greek small letter beta}'
            [''],
            ['LOA', self.length, 'm'],
            ['AHR', self.ahr, 'm, average hull roughness'],
            [''],
            ['---ATTITUDE---'],
            ['z_wl', self.z_wl, 'm, vertical distance of center of gravity to the calm water line'],
            ['tau', self.tau, 'deg, trim angle'],
            ['\u03B7\u2083', self.eta_3, 'deg, additional heave'],
            ['\u03B7\u2085', self.eta_5, 'deg, additional trim'],
            ['Transom draft', self.L_K*np.sin((self.tau+self.eta_5)*np.pi/180), 'm, draft of keel at transom'],
            [''],
            ['---PROPULSION---'],
            ['Thrust angle', self.epsilon, 'deg w.r.t. keel (CCW with body-fixed origin at 9 o\'clock)'],
            ['LCT', self.lT, 'm from stern, positive forward'],
            ['VCT', self.vT, 'm from keel, positive up'],
            [''],
            ['---FLAP---'],
            ['Chord', self.Lf, 'm'],
            ['Span/Beam', self.sigma, ''],
            ['Angle', self.delta, 'deg w.r.t. keel (CCW with body-fixed origin at 9 o\'clock)'],
            [''],
            ['---AIR DRAG---'],
            ['l_air', self.l_air, 'm, distance from stern to center of air pressure'],
            ['h_air', self.h_air, 'm, height from keel to top of square which bounds the air-drag-inducing shape'],
            ['b_air', self.b_air, 'm, transverse width of square which bounds the air-drag-inducing shape'],
            ['C_shape', self.C_shape, 'area coefficient for air-drag-inducing shape. C_shape = 1 means the air drag reference area is h_air*b_air'],
            ['C_D', self.C_D, 'air drag coefficient'],
            [''],
            ['---ENVIRONMENT---'],
            ['\u03C1', self.rho, 'kg/m\u00B3, water density'],
            ['\u03BD', self.nu, 'm\u00B2/s, water kinematic viscosity'],
            ['\u03C1_air', self.rho_air, 'kg/m\u00B3, air density'],
            ['g', self.g, 'm/s\u00B2, gravitational acceleration'],
            [''],
            ['---WETTED LENGTH OPTIONS---'],
            ['LC_type', self.LC_type, '(1 = Use Faltinsen 2010 wave rise approximation, 2 = Use Savitsky\'s \'64 approach, 3 = Use Savitsky\'s \'76 approach)'],
            ['zmax_type', self.zmax_type, '(1 = Uses 3rd order polynomial fit (faster, recommended), 2 = Use cubic interpolation)'],
            [''],
            ['---WETTED LENGTHS---'],
            ['L_K', self.L_K, 'm, keel wetted length'],
            ['L_C', self.L_C, 'm, chine wetted length'],
            ['\u03BB', self.lambda_W, 'mean wetted-length to beam ratio (L_K+L_C)/(2*beam)'],
            ['x_s', self.x_s, 'm, distance from keel/water-line intersection to start of wetted chine'],
            ['z_max', self.z_max, 'maximum presssure coordinate coefficient (z_max/Ut)'],
            [''],
            ['---FORCES [F_x (N, +aft), F_z (N, +up), M_cg (N*m, +pitch up)]---'],
            ['Hydrodynamic Force', self.hydrodynamic_force, ''],
            ['Skin Friction', self.skin_friction, ''],
            ['Air Resistance', self.air_resistance, ''],
            ['Flap Force', self.flap_force, ''],
            ['Net Force', self.net_force, ''],
            ['Resultant Thrust', self.thrust_force, ''],
            [''],
            ['---THURST & POWER---'],
            ['Thrust Magnitude', np.sqrt(self.thrust_force[0]**2+self.thrust_force[1]**2), 'N'],
            ['Effective Thrust', -self.thrust_force[0], 'N'],
            ['Eff. Power', -self.thrust_force[0]*self.speed/1000, 'kW'],
            ['Eff. Horsepower', -self.thrust_force[0]*self.speed/1000/0.7457, 'hp'],
            [''],
            ['---EOM MATRICES---'],
            ['Mass matrix, [kg, kg*m/rad; kg*m, kg*m\u00B2/rad]', self.mass_matrix, ''],
            ['Damping matrix, [kg/s, kg*m/(s*rad); kg*m/s, kg*m\u00B2/(s*rad)]', self.damping_matrix, ''],
            ['Restoring matrix, [N/m, N/rad; N, N*m/rad]', self.restoring_matrix, ''],
            [''],
            ['---PORPOISING---'],
            ['[[Eigenvalue check result, Est. pitch settling time (s)],\n [Savitsky chart result, Critical trim angle (deg)]]', np.array(self.porpoising), ''],
            [''],
            ['---BEHAVIOR IN WAVES---'],
            ['H_sig', self.H_sig, 'm, significant wave heigth'],
            ['R_AW', self.R_AW, 'N, added resistance in waves'],
            ['Average impact acceleration [n_cg, n_bow] (g\'s)', self.avgImpactAcc, ''],
        ]            
        
        cLens=[16,0,0] #Min spacing for columns
        for row in table:
            if len(row)==3:
                if row[1] is None:
                    print('{desc:<{cL0}} {val:<{cL1}} {unit:<{cL2}}'.format(desc=row[0], val=row[1], unit='None', cL0='', cL1=cLens[1], cL2=cLens[2]))
                elif isinstance(row[1], (list,np.ndarray)):
                    print(row[0]+' =')
                    with np.printoptions(formatter={'float': f'{{:.{sigFigs}g}}'.format}):
                        print(row[1])
                    print(row[2])
                else:
                    print('{desc:<{cL0}} {val:<{cL1}.{sNum}g} {unit:<{cL2}}'.format(desc=row[0], val=row[1], unit=row[2], cL0=cLens[0], cL1=cLens[1], cL2=cLens[2], sNum=sigFigs))
            else:
                print(row[0])
        
    def get_geo_lengths(self):
        """This function outputs the geometric lengths. 
        
        Adds/Updates the following parameters:
        --------------------------------------
        L_K (float): m, keel wetted length
        L_C (float): m, chine wetted length
        lambda_W (float): dimensionless, mean wetted-length to beam ratio (L_K+L_C)/(2*beam)
        x_s (float): m, distance from keel/water-line intersection to start of wetted chine
        z_max (float): dimensionless, maximum presssure coordinate coefficient (z_max/Ut)
        """
        b = self.beam
        lcg = self.lcg
        vcg = self.vcg
        z_wl = self.z_wl
        tau = self.tau
        beta = self.beta
        eta_3 = self.eta_3
        eta_5 = self.eta_5
        pi = np.pi
        wetted_lengths_type = self.LC_type
        z_max_type = self.zmax_type
        
        #Keel wetted length, Eq. 9.50 of Faltinsen 2010
        L_K = lcg + vcg / np.tan(pi/180*(tau + eta_5)) - (z_wl + eta_3) / np.sin(pi/180*(tau + eta_5))
        if L_K < 0:
            L_K = 0
        
        if wetted_lengths_type == 1:
            #z_max/Vt coefficient, Table 8.3 of Faltinsen 2010, page 303---------------
            beta_table = [4, 7.5, 10, 15, 20, 25, 30, 40]
            z_max_table = [0.5695, 0.5623, 0.5556, 0.5361, 0.5087, 0.4709, 0.4243, 0.2866]

            #Extrapolation warning
            if beta < beta_table[0] or beta > beta_table[-1]:
                warnings.warn('Deadrise ({0:.3f}) outside the interpolation range of 4-40 deg (Table 8.3 of Faltinsen 2010). Extrapolated values might be inaccurate.'.format(beta), stacklevel=2)

            if z_max_type == 1:
                z_max = np.polyval([-2.100644618790201e-006, -6.815747611588763e-005, -1.130563334939335e-003, 5.754510457848798e-001], beta)
            elif z_max_type == 2:
                z_max_func = interpolate.interp1d(beta_table, z_max_table, kind='cubic', fill_value='extrapolate') #Interpolation of the table
                z_max = z_max_func(beta)
            #--------------------------------------------------------------------------

            #Distance from keel/water-line intersection to start of wetted chine (Eq. 9.10 of Faltinsen)
            x_s = 0.5 * b * np.tan(pi/180*beta) / ((1 + z_max) * (pi/180)*(tau + eta_5))
            if x_s < 0:
                x_s = 0

            #Chine wetted length, Eq. 9.51 of Faltinsen 2010
            L_C = L_K - x_s
            if L_C < 0:
                L_C = 0
                x_s = L_K
                warnings.warn('Vessel operating with dry chines (L_C = 0).', stacklevel=2)

            #Mean wetted length-to-beam ratio
            lambda_W = (L_K + L_C) / (2 * b)
            
        elif wetted_lengths_type == 2:
            #Eq. 3 of Savitsky '64
            x_s = b/pi*np.tan(pi/180*beta)/np.tan(pi/180*(tau + eta_5))
            
            #Chine wetted length
            L_C = L_K - x_s
            if L_C < 0:
                L_C = 0
                x_s = L_K
                warnings.warn('Vessel operating with dry chines (L_C = 0).', stacklevel=2)
            
            #Mean wetted length-to-beam ratio
            lambda_W = (L_K + L_C)/(2*b)

            #z_max/Vt coefficient (E. 9.10 of Faltinsen 2010 rearranged)
            z_max = 0.5 * b * np.tan(pi/180*beta) / (x_s * (pi/180)*(tau + eta_5)) - 1
        
        elif wetted_lengths_type == 3:
            #Eq. 12 of Savitsky '76
            w = (0.57 + beta/1000)*(np.tan(pi/180*beta)/(2*np.tan(pi/180*(tau+eta_5)))-beta/167)

            lambda_K = L_K/b

            #Eq. 14 of Savitsky '76 
            lambda_C = (lambda_K-w)-0.2*np.exp(-(lambda_K-w)/0.3)
            if lambda_C < 0:
                lambda_C = 0
            L_C = lambda_C*b

            #Mean wetted length-to-beam ratio, Eq. 15 of Savitsky '76
            lambda_W = (lambda_K + lambda_C)/2+0.03

            x_s = L_K-L_C

            #z_max/Vt coefficient (E. 9.10 of Faltinsen 2010 rearranged)
            z_max = 0.5 * b * np.tan(pi/180*beta) / (x_s * (pi/180)*(tau + eta_5)) - 1
        
        if self.length is not None:
            if L_K > self.length:
                warnings.warn('The estimated wetted chine length ({0:.3f}) is larger than the vessel length ({1:.3f}).'.format(L_K, self.length), stacklevel=2)
        
        #Update values
        self.L_K = L_K
        self.L_C = L_C
        self.lambda_W = lambda_W
        self.x_s = x_s
        self.z_max = z_max
    
    def get_forces(self, runGeoLengths=True):
        """This function calls all the force functions to update the respective object attributes.
        
        PARAMETERS
        ----------
        runGeoLengths (boolean): Calculate the wetted lengths before calculating the forces.
        
        METHODS
        -------
        get_hydrodynamic_force()
        get_skin_friction()
        get_air_resistance()
        get_flap_force()
        sum_forces()
        
        Adds/updates the following parameters:
        --------------------------------------
        hydrodynamic_force (numpy float array): [F_x, F_z, M_cg] N, hydrodynamic force (x, y, rot directions in intertial coordinates)
        skin_friction (numpy float array): [F_x, F_z, M_cg] N, skin friction force
        air_resistance (numpy float array): [F_x, F_z, M_cg] N, air resistance force
        flap_force (numpy float array): [F_x, F_z, M_cg] N, flap resultant force
        thrust_force (numpy float array): [F_x, F_z, M_cg] N, thrust resultant force
        net_force (numpy float array): [F_x, F_z, M_cg] N, net force
        """
        if runGeoLengths:
            self.get_geo_lengths() #Calculated wetted lengths in get_forces()
        
        g = self.g 
        rho_air = self.rho_air
        C_D = self.C_D
        rho = self.rho
        nu = self.nu
        AHR = self.ahr
        W = self.weight
        epsilon = self.epsilon
        vT = self.vT
        lT = self.lT
        U = self.speed
        b = self.beam
        lcg = self.lcg
        vcg = self.vcg
        
        Lf = self.Lf
        sigma = self.sigma
        delta = self.delta
        beam = self.beam
        
        l_air = self.l_air
        h_air = self.h_air
        b_air = self.b_air
        C_shape = self.C_shape
        
        z_wl = self.z_wl
        tau = self.tau
        beta = self.beta
        eta_3 = self.eta_3
        eta_5 = self.eta_5
        
        L_K = self.L_K
        L_C = self.L_C
        lambda_W = self.lambda_W
        x_s = self.x_s
        z_max = self.z_max
        
        pi = np.pi
                
        def get_hydrodynamic_force():
            """This function follows Savitsky 1964 and Faltinsen 2005 in calculating the vessel's hydrodynamic forces and moment.
            """
            #Beam Froude number
            Fn_B = U/np.sqrt(g*b)
            
            #Warnings
            if Fn_B < 0.6 or Fn_B > 13:
                warnings.warn('Beam Froude number = {0:.3f}, outside of range of applicability (0.60 <= U/sqrt(g*b) <= 13.00) for planing lift equation. Results are extrapolations.'.format(Fn_B), stacklevel=2)
            if lambda_W > 4:
                warnings.warn('Mean wetted length-beam ratio = {0:.3f}, outside of range of applicability (lambda <= 4) for planing lift equation. Results are extrapolations.'.format(lambda_W), stacklevel=2)
            if tau < 2 or tau > 15:
                warnings.warn('Vessel trim = {0:.3f}, outside of range of applicability (2 deg <= tau <= 15 deg) for planing lift equation. Results are extrapolations.'.format(tau), stacklevel=2)

            #0-Deadrise lift coefficient
            C_L0 = (tau + eta_5)**1.1 * (0.012 * lambda_W**0.5 + 0.0055 * lambda_W**2.5 / Fn_B**2)

            #Lift coefficient with deadrise, C_Lbeta
            C_Lbeta = C_L0 - 0.0065 * beta * C_L0**0.6

            #Vertical force (lift)
            F_z = C_Lbeta * 0.5 * rho * U**2 * b**2

            #Horizontal force
            F_x = F_z*np.tan(pi/180*(tau + eta_5))

            #Lift's Normal force w.r.t. keel
            F_N = F_z / np.cos(pi/180*(tau + eta_5))

            #Longitudinal position of the center of pressure, l_p (Eq. 4.41, Doctors 1985)
            l_p = lambda_W * b * (0.75 - 1 / (5.21 * (Fn_B / lambda_W)**2 + 2.39)) #Limits for this is (0.60 < Fn_B < 13.0, lambda < 4.0)

            #Moment about CG (Axis consistent with Fig. 9.24 of Faltinsen (P. 366)
            M_cg = - F_N * (lcg - l_p)
            
            #Update values
            self.hydrodynamic_force = np.array([F_x, F_z, M_cg])
            
        def get_skin_friction():
            """This function outputs the frictional force of the vessel using ITTC 1957 and the Bowden and Davison 1974 roughness coefficient
            """
            #Surface area of the dry-chine region
            S1 = x_s * b / (2 * np.cos(pi/180*beta)) 
            if L_K < x_s:
                S1 = S1 * (L_K / x_s)**2

            #Surface area of the wetted-chine region
            S2 = b * L_C / np.cos(pi/180*beta) 

            #Total surface area
            S = S1 + S2 
            if S == 0: 
                F_x = 0
                F_z = 0
                M_cg = 0
            else:
                #Mean bottom fluid velocity, Savitsky 1964 - Hadler's empirical formula
                V_m = U * np.sqrt(1 - (0.012 * tau**1.1 * np.sqrt(lambda_W) - 0.0065 * beta * (0.012 * np.sqrt(lambda_W) * tau**1.1)**0.6) / (lambda_W * np.cos(tau * pi/180)))

                #Reynolds number (with bottom fluid velocity)
                Rn = V_m * lambda_W * b / nu

                #'Friction coefficient' ITTC 1957
                C_f = 0.075/(np.log10(Rn) - 2)**2

                #Additional 'friction coefficient' due to skin friction, Bowden and Davison (1974)
                deltaC_f = (44*((AHR/(lambda_W*b))**(1/3) - 10*Rn**(-1/3)) + 0.125)/10**3

                #Frictional force
                R_f = 0.5 * rho * (C_f + deltaC_f) * S * U**2

                #Geometric vertical distance from keel
                l_f = (b / 4 * np.tan(pi/180*beta) * S2 + b / 6 * np.tan(pi/180*beta) * S1) / (S1 + S2)

                #Horizontal force
                F_x = R_f * np.cos(pi/180*(tau + eta_5))

                #Vertical force
                F_z = - R_f * np.sin(pi/180*(tau + eta_5))

                #Moment about CG (Axis consistent with Fig. 9.24 of Faltinsen (P. 366)
                M_cg = R_f * (l_f - vcg)
                
            #Update values
            self.skin_friction = np.array([F_x, F_z, M_cg])
                
        def get_air_resistance():
            """This function estimates the air drag. It assumes a square shape projected area with a shape ceofficient.
            """
            if C_shape == 0 or b_air == 0:
                self.air_resistance = np.array([0, 0, 0])
                return

            #Vertical distance from calm water line to keel at LOA
            a_dist = np.sin(pi/180*(tau + eta_5))*(l_air-L_K)
            
            #Vertical distance from keel to horizontal line level with ACV's height
            b_dist = np.cos(pi/180*(tau + eta_5))*h_air
            
            #Vertical distance from CG to center of square (moment arm, positive is CG above)
            momArm = z_wl - (a_dist + b_dist)/2
            
            #Square projected area
            Area = (a_dist+b_dist)*b_air*C_shape
            if Area < 0:
                Area = 0
                
            #Horizontal force (Positive aft)
            F_x = 0.5*rho_air*C_D*Area*U**2
            
            #Vertical force (Positive up) 
            F_z = 0
            
            #Moment (positve CCW)
            M_cg = -F_x*momArm

            #Update values
            self.air_resistance = np.array([F_x, F_x, M_cg])
        
        def get_flap_force():
            """This function outputs the flap forces w.r.t. global coordinates (Savitsky & Brown 1976). Horz: Positive Aft, Vert: Positive Up, Moment: Positive CCW
            """
            if Lf == 0:
                self.flap_force = np.array([0, 0, 0])
                return
            
            #Warnings
            if Lf > 0.10*(L_K + L_C)/2 or Lf < 0:
                warnings.warn('Flap chord = {0:.3f} outside of bounds (0-10% of mean wetted length) for flap forces estimates with Savitsky & Brown 1976'.format(Lf), stacklevel=2)
            if delta < 0 or delta > 15:
                warnings.warn('Flap deflection angle = {0:.3f} out of bounds (0-15 deg) for flap forces estimates with Savitsky & Brown 1976'.format(delta), stacklevel=2)
            Fn_B = U/np.sqrt(g*b)
            if Fn_B < 2 or Fn_B > 7:
                warnings.warn('Beam-based Froude number Fn_B = {0:.3f} out of bounds (2-7) for flap forces estimates with Savitsky & Brown 1976'.format(Fn_B), stacklevel=2)
            
            F_z = 0.046*(Lf*3.28084)*delta*sigma*(b*3.28084)*(rho/515.379)/2*(U*3.28084)**2*4.44822

            F_x = 0.0052*F_z*(tau+eta_5+delta)

            l_flap = 0.6*b+Lf*(1-sigma)

            M_cg = -F_z*(lcg-l_flap)
            
            #Update values
            self.flap_force = np.array([F_x, F_z, M_cg])
        
        def sum_forces():
            """This function gets the sum of forces and moments, and consequently the required net thrust. The coordinates are positive aft, positive up, and positive counterclockwise.
            """
            #Call all force functions-------
            get_hydrodynamic_force()
            get_skin_friction()
            get_air_resistance()
            get_flap_force()
            #-------------------------------
            
            forcesMatrix = np.column_stack((self.gravity_force, self.hydrodynamic_force, self.skin_friction, self.air_resistance, self.flap_force)) #Forces and moments
            F_sum = np.sum(forcesMatrix, axis=1) #F[0] is x-dir, F[1] is z-dir, and F[2] is moment

            #Required thrust and resultant forces
            T = F_sum[0]/np.cos(pi/180*(epsilon+tau+eta_5)); #Magnitude
            T_z = T*np.sin(pi/180*(epsilon+tau+eta_5)); #Vertical
            T_cg = T*np.cos(pi/180*epsilon)*(vcg - vT) - T*np.sin(pi/180*epsilon)*(lcg - lT); #Moment about cg
            
            #Update resultant thurst values
            self.thrust_force = np.array([-F_sum[0], T_z, T_cg])
            
            #Include resultant thrust forces in sum
            F_sum[1] = F_sum[1]+T_z
            F_sum[2] = F_sum[2]+T_cg
            
            #Update values
            self.net_force = F_sum
            
        #Call functions
        sum_forces()
        
    def get_steady_trim(self, x0=[0, 3], tauLims=[0.5, 35], tolF=10**-6, maxiter=50):
        """This function finds and sets the equilibrium point when the vessel is steadily running in calm water.
        
        PARAMETERS
        ----------
        x0 (float array): Initial guess for equilibirum point [z_wl (m), tau (deg)] (default of [0, 3])
        tauLims (float array): Limits for equilibrium trim search (default of [0.5, 35])
        tolF (float): Tolerance for convergence to zero (default of 10**-6)
        maxiter (float): Maximum iterations (default of 50)
        
        Updates the following parameters:
        --------------------------------------
        z_wl (float): m, vertical distance of center of gravity to the calm water line (default of 0)
        tau (float): deg, trim angle (default of 5)
        """
        def _boatForces(x):
            self.z_wl = x[0]
            self.tau = x[1]
            self.get_forces()
            return self.net_force[1:3]

        def _boatForcesPrime(x):
            return ndmath.complexGrad(_boatForces, x)

        def _L_K(x):
            self.z_wl = x[0]
            self.tau = x[1]
            self.get_geo_lengths()
            return [-self.L_K]
        
        xlims = np.array([[-np.Inf, np.Inf], tauLims])
        warnings.filterwarnings("ignore", category=UserWarning)
        _ = ndmath.nDimNewton(_boatForces, x0, _boatForcesPrime, tolF, maxiter, xlims, hehcon=_L_K)
        warnings.filterwarnings("default", category=UserWarning)
        
    def get_eom_matrices(self, runGeoLengths=True):
        """This function returns the mass, damping, and stiffness matrices following Faltinsen 2010.
        
        PARAMETERS
        ----------
        runGeoLengths (boolean): Calculate the wetted lengths before calculating the EOM matrices (default of True)
        
        METHODS
        -------
        get_mass_matrix()
        get_damping_matrix()
        get_restoring_matrix(diffType=1, step=10**-6.6)
        
        Adds/updates the following parameters:
        --------------------------------------
        mass_matrix (numpy float 2D array): [[A_33 (kg), A_35 (kg*m/rad)], [A_53 (kg*m), A_55 (kg*m^2/rad)]]
        damping_matrix (numpy float 2D array): [[B_33 (kg/s), B_35 (kg*m/(s*rad))], [B_53 (kg*m/s), B_55 (kg*m**2/(s*rad))]]
        restoring_matrix (numpy float 2D array): [[C_33 (N/m), C_35 (N/rad)], [C_53 (N), C_55 (N*m/rad)]]
        """
        if runGeoLengths:
            self.get_geo_lengths() #Calculated wetted lengths in get_eom_matrices()
        
        W = self.weight
        U = self.speed
        rho = self.rho
        b = self.beam
        lcg = self.lcg
        tau = self.tau
        beta = self.beta
        g = self.g
        r_g = self.r_g
        
        eta_5 = self.eta_5
        
        L_K = self.L_K
        L_C = self.L_C
        lambda_W = self.lambda_W
        x_s = self.x_s
        z_max = self.z_max
        
        pi = np.pi
        
        def get_mass_matrix():
            """This function returns the added mass coefficients following Sec. 9.4.1 of Faltinsen 2010, including weight and moment of inertia
            """
            
            #Distance of CG from keel-WL intersection
            x_G = L_K - lcg

            #K constant (Eq. 9.63 of Faltinsen 2010)
            K = (pi / np.sin(pi/180*beta) * gamma(1.5 - beta/180) / (gamma(1 - beta/180)**2 * gamma(0.5 + beta/180)) - 1) / np.tan(pi/180*beta)

            kappa = (1 + z_max) * (pi/180)*(tau + eta_5) #User defined constant

            #Based on Faltinsen's
            A1_33 = rho * kappa**2 * K * x_s**3 / 3
            A1_35 = A1_33 * (x_G - x_s * 3/4)
            A1_53 = A1_35
            A1_55 = A1_33 * (x_G**2 - 3/2 * x_G * x_s + 3/5 * x_s**2)

            #Contribution from wet-chine region
            if L_C > 0:
                C_1 = 2 * np.tan(pi/180*beta)**2 / pi * K

                A2_33 = (rho * b**3) * C_1 * pi / 8 * L_C / b
                A2_35 = (rho * b**4) * (- C_1 * pi / 16 * ((L_K / b)**2 - (x_s / b)**2) + x_G / b * A2_33 / (rho * b**3))
                A2_53 = A2_35
                A2_55 = (rho * b**5) * (C_1 * pi / 24 * ((L_K / b)**3 - (x_s / b)**3) - C_1 / 8 * pi * (x_G / b) * ((L_K / b)**2 - (x_s / b)**2) + (x_G / b)**2 * A2_33 / (rho * b**3))
            else:
                A2_33 = 0
                A2_35 = 0
                A2_53 = 0
                A2_55 = 0

            #Total added mass & update values
            A_33 = A1_33 + A2_33 + W/g # kg, A_33
            A_35 = A1_35 + A2_35 # kg*m/rad, A_35
            A_53 = A1_53 + A2_53 # kg*m, A_53
            A_55 = A1_55 + A2_55 + W/g*r_g**2 # kg*m^2/rad, A_55
            self.mass_matrix = np.array([[A_33, A_35], [A_53, A_55]])
            
        def get_damping_matrix():
            """This function returns the damping coefficients following Sec. 9.4.1 of Faltinsen 2010
            """
            #Heave-heave added mass (need to substract W/g since it was added)
            A_33 = self.mass_matrix[0,0] - W/g

            if L_C > 0:
                d = 0.5 * b * np.tan(pi/180*beta)
            else:
                d = (1 + z_max) * (pi/180)*(tau + eta_5) * L_K

            #K constant (Eq. 9.63 of Faltinsen 2010, P. 369)
            K = (pi / np.sin(pi/180*beta) * gamma(1.5 - beta/180) / (gamma(1 - beta/180)**2 * gamma(0.5 + beta/180)) - 1) / np.tan(pi/180*beta)

            #2D Added mass coefficient in heave
            a_33 = rho * d**2 * K

            #Infinite Fn lift coefficient
            C_L0 = (tau + eta_5)**1.1 * 0.012 * lambda_W**0.5

            #Derivative w.r.t. tau (rad) of inf. Fn C_L0
            dC_L0 = (180 / pi)**1.1 * 0.0132 * (pi/180*(tau + eta_5))**0.1 * lambda_W**0.5

            #Derivative w.r.t. tau (rad) of inf. Fn C_Lbeta
            dC_Lbeta = dC_L0 * (1 - 0.0039 * beta * C_L0**-0.4)

            #Damping coefficients & update values
            B_33 = rho / 2 * U * b**2 * dC_Lbeta # kg/s, B_33, Savitsky based
            B_35 = - U * (A_33 + lcg * a_33) # kg*m/(s*rad), B_35, Infinite frequency based
            B_53 = B_33 * (0.75 * lambda_W * b - lcg) # kg*m/s, B_53, Savitsky based
            B_55 = U * lcg**2 * a_33 # kg*m**2/(s*rad), B_55, Infinite frequency based
            self.damping_matrix = np.array([[B_33, B_35], [B_53, B_55]])
            
        def get_restoring_matrix(diffType=1, step=10**-6.6):
            """This function returns the restoring coefficients following the approach in Sec. 9.4.1 of Faltinsen 2010
            
            PARAMETERS
            ----------
            diffType (int): 1 (default, recommended) = Complex step method, 2 = Foward step difference
            step (float): Step size if using diffType == 2 (default of 10**-6)
            """
            def func(eta):
                self.eta_3 = eta[0] 
                self.eta_5 = eta[1]
                self.get_forces()
                return self.net_force[1:3]
            
            temp_eta_3 = self.eta_3
            temp_eta_5 = self.eta_5
            
            if diffType == 1:
                C_full = -ndmath.complexGrad(func, [temp_eta_3, temp_eta_5])
            elif diffType == 2:
                C_full = -ndmath.finiteGrad(func, [temp_eta_3, temp_eta_5], 10**-6.6)

            #Reset values
            self.eta_3 = temp_eta_3
            self.eta_5 = temp_eta_5
            self.get_forces()
            
            #Conversion deg to rad (degree in denominator)
            C_full[0,1] = C_full[0,1] / (pi/180) # N/rad, C_35
            C_full[1,1] = C_full[1,1] / (pi/180) # N*m/rad, C_55
            
            #Update values
            self.restoring_matrix = C_full
        
        #Call functions
        get_mass_matrix()
        get_damping_matrix()
        get_restoring_matrix()
    
    def check_porpoising(self, stepEstimateType=1):
        """This function checks for porpoising.
        
        PARAMETERS
        ----------
        stepEstimateType (integer): pitch step response settling time estimate type, 1 (default) = -3/np.real(eigVals[0])], 2 = time-domain simulation estimate
        
        Adds/updates the following parameters
        -------------------------------------
        porpoising (list): [[eigenvalue result, est. pitch settling time (s)], [Savitsky chart result, critical trim angle (deg)]]
        """
        #Eigenvalue analysis
        try:
            self.mass_matrix
        except AttributeError:
            warnings.warn('No Equation Of Motion (EOM) matrices found. Running get_eom_matrices().', stacklevel=2)
            self.get_eom_matrices()
            
        M = self.mass_matrix
        C = self.damping_matrix
        K = self.restoring_matrix
        
        nDim = len(M)
        A_ss = np.concatenate((np.concatenate((np.zeros((nDim,nDim)), np.identity(nDim)), axis=1), np.concatenate((-np.linalg.solve(M,K), -np.linalg.solve(M,C)), axis=1))) #State space reprecentation 
        
        eigVals = np.linalg.eigvals(A_ss)

        eig_porpoise = any(eigVal >= 0 for eigVal in eigVals)
        
        if stepEstimateType == 1:
            settling_time = -3/np.real(eigVals[0])
        elif stepEstimateType == 2: 
            B_ss = np.array([[1],[0],[0],[0]]) #Pitch only
            C_ss = np.array([[1,0,0,0]]) #Pitch only
            D_ss = np.array([[0]])
            
            system = (A_ss,B_ss,C_ss,D_ss)
            t, y = signal.step(system)
            settling_time = (t[next(len(y)-i for i in range(2,len(y)-1) if abs(y[-i]/y[-1])>1.02)]-t[0])
        
        #Savitsky '64 chart method
        C_L = self.weight/(1/2*self.rho*self.speed**2*self.beam**2)
        x = np.sqrt(C_L/2)

        #Warnings
        if x > 0.3 or x < 0.13:
            warnings.warn('Lift Coefficient = {0:.3f} outside of bounds (0.0338-0.18) for porpoising estimates with Savitsky 1964. Results are extrapolations.'.format(C_L), stacklevel=2)
        if self.beta > 20:
            warnings.warn('Deadrise = {0:.3f} outside of bounds (0-20 deg) for porpoising estimates with Savitsky 1964. Results are extrapolations.'.format(self.beta), stacklevel=2)

        tau_crit_0 = -376.37*x**3 + 329.74*x**2 - 38.485*x + 1.3415
        tau_crit_10 = -356.05*x**3 + 314.36*x**2 - 41.674*x + 3.5786
        tau_crit_20 = -254.51*x**3 + 239.65*x**2 - 23.936*x + 3.0195

        tau_crit_func = interpolate.interp1d([0, 10, 20], [tau_crit_0, tau_crit_10, tau_crit_20], kind='quadratic', fill_value='extrapolate')
        tau_crit = tau_crit_func(self.beta)

        if self.tau > tau_crit:
            chart_porpoise = True
        else:
            chart_porpoise = False
        
        #Update values
        self.porpoising = [[eig_porpoise, settling_time], [chart_porpoise, float(tau_crit)]]
        
    def get_seaway_behavior(self):
        """This function calculates the seaway behavior as stated in Savitsky & Brown '76.
        
        Adds/updates the following parameters
        -------------------------------------
        avg_impact_acc (numpy float array): [n_cg, n_bow] g's, average impact acceleration at center of gravity and bow
        R_AW (float): N, added resistance in waves
        """
        if self.H_sig is None:
            self.H_sig = self.beam*0.5 #Arbitrary wave height if no user-defined wave height
            warnings.warn('Significant wave height has not been specified. Using beam*0.5 = {0:.3f} m.'.format(self.H_sig), stacklevel=2)
        if self.length is None:
            self.length = self.beam*3
            warnings.warn('Vessel length has not been specified. Using beam*3 = {0:.3f} m.'.format(self.length), stacklevel=2)
        H_sig = self.H_sig
        
        W = self.weight
        beta = self.beta
        tau = self.tau
        
        pi = np.pi
        
        Delta_LT = W/9964 #Displacement in long tons
        Delta = Delta_LT*2240 #Displacement in lbf
        L = self.length*3.281 #Length in ft
        b = self.beam*3.281 #Beam in ft
        V_K = self.speed*1.944 #Speed in knots
        H_sig = H_sig*3.281 #Significant wave height in ft
        
        w = self.rho*self.g/(4.448*35.315) #Specific weight in lbf/ft^3
        
        C_Delta = Delta/(w*b**3) #Static beam-loading coefficient
               
        #Check that variables are inside range of applicability (P. 395 of Savitsky & Brown '76)
        P1 = Delta_LT/(0.01*L)**3
        P2 = L/b
        P5 = H_sig/b
        P6 = V_K/np.sqrt(L)
        if P1 < 100 or P1 > 250:
            warnings.warn('Vessel displacement coefficient = {0:.3f}, outside of range of applicability (100 <= Delta_LT/(0.01*L)^3 <= 250, with units LT/ft^3). Results are extrapolations.'.format(P1), stacklevel=2)
        if P2 < 3 or P2 > 5:
            warnings.warn('Vessel length/beam = {0:.3f}, outside of range of applicability (3 <= L/b <= 5). Results are extrapolations.'.format(P2), stacklevel=2)
        if tau < 3 or tau > 7:
            warnings.warn('Vessel trim = {0:.3f}, outside of range of applicability (3 deg <= tau <= 7 deg). Results are extrapolations.'.format(tau), stacklevel=2)
        if beta < 10 or beta > 30:
            warnings.warn('Vessel deadrise = {0:.3f}, outside of range of applicability (10 deg <= beta <= 30 deg). Results are extrapolations.'.format(beta), stacklevel=2)
        if P5 < 0.2 or P5 > 0.7:
            warnings.warn('Significant wave height / beam = {0:.3f}, outside of range of applicability (0.2 <= H_sig/b <= 0.7). Results are extrapolations.'.format(P5), stacklevel=2)
        if P6 < 2 or P6 > 6:
            warnings.warn('Speed coefficient = {0:.3f}, outside of range of applicability (2 <= V_K/sqrt(L) <= 6, with units knots/ft^0.5). Results are extrapolations.'.format(P6), stacklevel=2)
            
        R_AW_2 = (w*b**3)*66*10**-6*(H_sig/b+0.5)*(L/b)**3/C_Delta+0.0043*(tau-4) #Added resistance at V_K/sqrt(L) = 2
        R_AW_4 = (Delta)*(0.3*H_sig/b)/(1+2*H_sig/b)*(1.76-tau/6-2*np.tan(beta*pi/180)**3) #V_K/sqrt(L) = 4
        R_AW_6 = (w*b**3)*(0.158*H_sig/b)/(1+(H_sig/b)*(0.12*beta-21*C_Delta*(5.6-L/b)+7.5*(6-L/b))) #V_K/sqrt(L) = 6
        R_AWs = np.array([R_AW_2, R_AW_4, R_AW_6])*4.448 #lbf to N conversion
        
        R_AWs_interp = interpolate.interp1d([2,4,6], R_AWs, kind='quadratic', fill_value='extrapolate')
        R_AW = R_AWs_interp([V_K/np.sqrt(L)])
        
        n_cg = 0.0104*(H_sig/b+0.084)*tau/4*(5/3-beta/30)*(V_K/np.sqrt(L))**2*L/b/C_Delta #g, at CG
        n_bow = n_cg*(1+3.8*(L/b-2.25)/(V_K/np.sqrt(L))) #g, at bow
        avg_impact_acc = np.array([n_cg, n_bow])
        
        #Update values
        self.avgImpactAcc = avg_impact_acc
        self.R_AW = R_AW[0]