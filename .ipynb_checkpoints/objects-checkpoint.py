import numpy as np, math as m

###CLASS
class solar_system:

    def __init__(self):

        self.ua = 150E9
        self.mO = 1.989E30
        self.G = 6.674184E-11
        self.x = 0
        self.y = 0
        self.corps = []


class corps:

    def __init__(self, type, a0, e, mass, nom, color, universe):

        self.name = nom
        self.mass = mass
        self.radius = a0
        self.exc = e
        self.type = type
        self.color = color

        self.a = self.radius * universe.ua
        if self.type == 1:
            self.period = 2*m.pi*m.sqrt(self.a**3/(universe.G*universe.mO))

        self.perihelie = self.a*(1-self.exc)
        self.aphelie = (2*self.a)-self.perihelie
        self.b = self.a*m.sqrt(1-self.exc**2)

        self.xi = self.perihelie
        self.yi = 0
        self.vxi = 0

        self.x = self.perihelie
        self.y = 0
        if self.type == 1:
            self.vx = 0
            self.vyi = m.sqrt(universe.G*universe.mO*(1+self.exc)/(self.a*(1-self.exc)))
            self.vy = self.vyi

        self.forces = []
        self.resultante = np.array([0., 0.]) #[Fx, Fy]

    def vecteur_directeur(self, ssystem):
        '''Vecteur position selon le soleil'''
        return np.array([self.x - ssystem.x, self.y - ssystem.y])


    def distance_corps(self, planete):

        return m.sqrt((self.x - planete.x)**2 + (self.y - planete.y)**2)

    def direction_corps(self, planete, ssystem):
        '''Vecteur direceteur entre les 2 corps'''
        v1, v2 = self.vecteur_directeur(ssystem), planete.vecteur_directeur(ssystem)
        return (v2 - v1)/self.distance_corps(planete)

    def refresh_cin(self):

        self.x = self.xi
        self.y = self.yi
        self.vx = self.vxi
        self.vy = self.vyi


class force:

    def __init__(self, beta, corps1, corps2, ssystem):

        self.k = ssystem.G * ssystem.mO
        self.G = ssystem.G
        self.beta = beta
        self.corps1 = corps1
        self.corps2 = corps2
        self.fx = 0
        self.fy = 0
        if self.corps1.name != "ref":
            self.dir = corps1.direction_corps(corps2, ssystem)


    def exp_force(self):

        return -self.corps1.mass*self.corps2.mass*self.G/((self.corps1.distance_corps(self.corps2))**self.beta)





###PLANETES
def init_systeme():

    systeme_solaire = solar_system()
    systeme_solaire_circ = solar_system()

    liste_planetes = [
                        corps(1, 0.39, 0.206, 2.4E23, "Mercury", "grey", systeme_solaire),
                        corps(1, 0.72, 0.007, 4.9E24, "Venus", "mediumseagreen", systeme_solaire),
                        corps(1, 1, 0.017, 6E24, "Earth", "royalblue", systeme_solaire),
                        corps(1, 1.52, 0.093, 6.6E23, "Mars", "orangered", systeme_solaire),
                        corps(1, 5.2, 0.048, 1.9E27, "Jupiter", "tan", systeme_solaire),
                        corps(1, 9.54, 0.056, 5.7E26, "Saturn", "gold", systeme_solaire),
                        corps(1, 19.19, 0.046, 8.8E25, "Uranus", "turquoise", systeme_solaire),
                        corps(1, 30.06, 0.01, 1E26, "Neptune", "mediumblue", systeme_solaire),
                        corps(1, 39.26, 0.248, 1.3E22, "Pluto", "sienna", systeme_solaire),
                    ]

    liste_planetes_circ = [
                        corps(1, 0.39, 0, 2.4E23, "Mercury", "grey", systeme_solaire),
                        corps(1, 0.72, 0, 4.9E24, "Venus", "mediumseagreen", systeme_solaire),
                        corps(1, 1, 0, 6E24, "Earth", "royalblue", systeme_solaire),
                        corps(1, 1.52, 0, 6.6E23, "Mars", "orangered", systeme_solaire),
                        corps(1, 5.2, 0, 1.9E27, "Jupiter", "tan", systeme_solaire),
                        corps(1, 9.54, 0, 5.7E26, "Saturn", "gold", systeme_solaire),
                        corps(1, 19.19, 0, 8.8E25, "Uranus", "turquoise", systeme_solaire),
                        corps(1, 30.06, 0, 1E26, "Neptune", "mediumblue", systeme_solaire),
                        corps(1, 39.26, 0, 1.3E22, "Pluto", "sienna", systeme_solaire),
                    ]

    soleil = corps(0, 0, 0, 1.989E30, "Sun", "orange", systeme_solaire)

    systeme_solaire.corps = [soleil, liste_planetes]
    systeme_solaire_circ.corps = [soleil, liste_planetes_circ]


    return (systeme_solaire, systeme_solaire_circ)