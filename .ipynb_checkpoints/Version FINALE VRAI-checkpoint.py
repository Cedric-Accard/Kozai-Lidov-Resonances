import os
import math as m
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from scipy.integrate import odeint
import numpy as np


os.chdir("C:\\Users\\proprietaire\\Desktop")
#os.chdir("C:\\Users\\elodie\\Desktop\\PMI")

from objects import *


### INIT
systeme_solaire, systeme_solaire_circ = init_systeme()
liste_planetes = systeme_solaire.corps[1]
liste_planetes_circ = systeme_solaire_circ.corps[1]
soleil = systeme_solaire.corps[0]

astre_ref = corps(0, 0, 0 ,0, "ref", "r", systeme_solaire)
fref = force(2, astre_ref, astre_ref, systeme_solaire)
k = fref.k

### REVOLUTION
def revolution(id_planet, beta, nb_periodes, dtau):

    '''
    Considère une planète, un beta et un temps d'intégration/pas et renvoie la liste des positions et des vitesses
    '''

    planet = liste_planetes[id_planet]
    U = (planet.x, planet.y, planet.vx, planet.vy)

    fg = force(beta, soleil, planet, systeme_solaire)
    k = fg.k

    #TRAJECTOIRE
    tmin = 0
    tps = nb_periodes * planet.period
    #tps = 5E13
    nb_points = int(tps/dtau)
    temps = np.linspace(tmin, tps, nb_points)

    def eq_mouvement(U,temps):

        x,y,dx,dy = U
        #cos_phi = (x+a*e)/a
        #r = a*(1-e*cos_phi) FAUX SI LA TRAJ N EST PAS ELLIPTIQUE (BETA != 2)
        r = m.sqrt(x**2 + y**2)

        return [dx,dy,-k*x/(r**(beta+1)),-k*y/(r**(beta+1))]


    X,Y,dX,dY = odeint(eq_mouvement, U, temps).T

    return (planet, temps, X, Y, dX, dY)


def print_trajectoire(planet, beta, X, Y):
    '''
    Pour une planète(classe), un beta et les positions en X et Y donnés, cette fonction affiche la trajectoire correspondante sur un plot
    '''
    a = planet.a
    e = planet.exc
    b = planet.b
    d = planet.perihelie
    q = planet.aphelie

    
    plt.figure()
    plt.plot(X,Y, '.')
    plt.plot(0,0,'o', -a*e,0,'x')
    plt.plot(-q,0,'.', d,0,'.')
    plt.title(planet.name + " beta=" + str(beta))
    plt.grid()
    plt.show()

def print_trajectoires(liste, beta, liste_cin):
    '''idem que la fonction précédente mais avec une liste de planetes, liste de X et de Y'''

    for p in range(len(liste)):
        planete = liste[p]
        print_trajectoire(planete, beta, liste_cin[p][0], liste_cin[p][1])
        
def print_trajectoire_dtauodeint(planet, beta, X, Y, texte):
    '''
    Même fonction que print_trajectoire() sauf que ici on peut customiser le titre pour que les graphes tracés soit adpatés à la méthode utilisé (nous odeint)
    '''
    a = planet.a
    e = planet.exc
    b = planet.b
    d = planet.perihelie
    q = planet.aphelie

    
    plt.figure()
    plt.plot(X,Y, '.')
    plt.plot(0,0,'o', -a*e,0,'x')
    plt.plot(-q,0,'.', d,0,'.')
    plt.title(texte)
    plt.grid()
    plt.show()




### SIMULATION ET ETUDE ENERGETIQUE POUR UNE PLANETE EN 2D
def simulation(id_planete, beta, nb_periodes, dtau):
    '''Pour une planète, fait le tracé de la trajectoire et l'étude énergétique
    '''
    planet, temps, X, Y, dX, dY = revolution(id_planete, beta, nb_periodes, dtau)
    print_trajectoire(planet, beta, X, Y)

    E, R, EC, EP, temps, Ecarts, Eth = energie(planet, beta, X, Y, dX, dY, temps)
    etude_energie(planet, E, R, EC, EP, Ecarts, Eth, temps)

### SIMULATION PLUSIEURS CORPS

#FONCTIONS AUXILIAIRES
def init_system(liste_corps):
    '''
    Initialise les forces (sous forme de liste) que subit chaque planète
    '''
    nb_corps = len(liste_corps)
    for p in liste_corps:
        p.forces = [force(2, p, astre_ref, systeme_solaire) for i in range(nb_corps)] #ref liste forces

def force_resultante(planete):
    '''
    Fonction qui pour chaque planète renvoie la force résultante des toutes les forces auxquelles est soumis une planète sous la forme d'un array 2D numpy """
    '''
    F = np.array([0., 0.])
    for f in planete.forces:

        if f.corps2.name != "ref":

            f.fx = f.exp_force() * np.dot(f.dir, np.array([1, 0]))
            f.fy = f.exp_force() * np.dot(f.dir, np.array([0, 1]))

            F[0] += f.fx
            F[1] += f.fy

    return F


## ANIMATIONS 2D

#TEMPS REEL
def systeme_temps_reel_2D(liste_corps, beta, tmax, dt, freq):
    '''
    Calcul et affiche,grâce à Euler, pas par pas, la trajectoire des planètes en 2D
    '''
    nb_planetes = len(liste_corps)
    fref = force(2, astre_ref, astre_ref, systeme_solaire)
    k = fref.k
    temps = 0
    max = liste_corps[-1].a
    fig = plt.figure()
    plt.ion()

    while temps < tmax:
        
        for p in range(nb_planetes):

            planete = liste_corps[p]
            fs = force(beta, soleil, planete, systeme_solaire)
            planete.forces[0] = fs

            for k in range(p+1, nb_planetes, 1):

                planete2 = liste_corps[k]
                f12 = force(beta, planete, planete2, systeme_solaire)
                f21 = force(beta, planete2, planete, systeme_solaire)

                planete.forces[k] = f21
                planete2.forces[p+1] = f12

        for p in liste_corps:
            p.resultante = force_resultante(p)
        
        for p in liste_corps:
            p.vx = p.vx + (p.resultante[0]/p.mass)*dt
            p.vy = p.vy + (p.resultante[1]/p.mass)*dt

            p.x = p.x + p.vx*dt
            p.y = p.y + p.vy*dt

            plt.plot([p.x], [p.y], '.', color = p.color)
        
        plt.plot([soleil.x], [soleil.y], 'o', color = soleil.color)
        plt.axis([-1.1*max,1.1*max,-1.1*max,1.1*max])
        
        plt.pause(1/freq)
        plt.show()
        plt.clf()

        temps+=dt
        

def simulation_temps_reel_2D(n, beta, temps, dt, freq):
    '''
    Fait une simulation pas par pas en utilisant init_system(liste) et
    systeme_temps_reel(liste, beta, temps, dt, freq), pour un groupe de planète, soit les 4 planètes intérieures, soit les 4 planètes extérieures. (selon valeur de n)
    '''
    if n == 0:
        liste = liste_planetes[0:4]
    elif n == 1:
        liste = liste_planetes[4:10]
    else:
        liste = liste_planetes
        
    for p in liste:
        p.refresh_cin()
        
    init_system(liste)
    systeme_temps_reel_2D(liste, beta, temps, dt, freq)



#PRECALCUL
def systeme_precalc_fast(liste_corps, beta, tmax, dtau, tmin=0):
    '''
    Sert à intégrer les equations différentielles pour chaque planète pas par pas. (Méthode d'Euler) (Car à chaque pas la force résultante change)
    Renvoie une liste contenant des sous listes contenant les positions X et Y pour chaque planète.
    '''
    nb_planetes = len(liste_corps)
    fref = force(2, astre_ref, astre_ref, systeme_solaire)
    k = fref.k

    

    tps = tmax - tmin
    nb_points = int(tps/dtau)
    temps = np.linspace(tmin, tps, nb_points)
        
    cinetique_all = [ [[], []] for p in liste_corps ]

    def acc(planete):

                Fx = planete.resultante[0]
                Fy = planete.resultante[1]

                return np.array([Fx/planete.mass, Fy/planete.mass])

    for t in temps:
        for p in range(nb_planetes):

            planete = liste_corps[p]
            fs = force(beta, soleil, planete, systeme_solaire)
            planete.forces[0] = fs #force liée au soleil

            for k in range(p+1, nb_planetes, 1):

                planete2 = liste_corps[k]
                f12 = force(beta, planete, planete2, systeme_solaire)
                f21 = force(beta, planete2, planete, systeme_solaire)

                planete.forces[k] = f21
                planete2.forces[p+1] = f12


            for m in range(nb_planetes):

                p = liste_corps[m]
                p.resultante = force_resultante(p)

                A = acc(p)

                p.vx += A[0]*dtau
                p.vy += A[1]*dtau
                p.x += p.vx*dtau
                p.y += p.vy*dtau
                

                cinetique_all[m][0].append(p.x)
                cinetique_all[m][1].append(p.y)

    return cinetique_all
            

def systeme_precalc_RK4(liste_corps, beta, tmax, dtau, tmin=0):
    '''
    Sert à intégrer les equations différentielles pour chaque planète pas par pas. (Méthode RK4)) (Car à chaque pas la force résultante change)
    Renvoie une liste contenant des sous listes contenant les positions X et Y pour chaque planète.
    '''
    nb_planetes = len(liste_corps)
    fref = force(2, astre_ref, astre_ref, systeme_solaire)
    k = fref.k
    tps = tmax - tmin
    nb_points = int(tps/dtau)
    temps = np.linspace(tmin, tps, nb_points)
        
    cinetique_all = [ [[], []] for p in liste_corps ]

    def eq_mouvement(U, planete):

                x, y, dx, dy = U
                Fx = planete.resultante[0]
                Fy = planete.resultante[1]

                return np.array([dx, dy, Fx/planete.mass, Fy/planete.mass])

    for t in temps:
        for p in range(nb_planetes):

            planete = liste_corps[p]
            fs = force(beta, soleil, planete, systeme_solaire)
            planete.forces[0] = fs #force liée au soleil

            for k in range(p+1, nb_planetes, 1):

                planete2 = liste_corps[k]
                f12 = force(beta, planete, planete2, systeme_solaire)
                f21 = force(beta, planete2, planete, systeme_solaire)

                planete.forces[k] = f21
                planete2.forces[p+1] = f12

            for m in range(nb_planetes):

                p = liste_corps[m]
                p.resultante = force_resultante(p)


                U = np.array([p.x, p.y, p.vx, p.vy])
                k1 = dtau * eq_mouvement(U, p)
                k2 = dtau * eq_mouvement(U+(k1/2), p)
                k3 = dtau * eq_mouvement(U+(k2/2), p)
                k4 = dtau * eq_mouvement(U+k3, p)
                k = (k1+2*k2+2*k3+k4)/6
                Un = U + k

                p.x = Un[0]
                p.y = Un[1]
                p.vx = Un[2]
                p.vy = Un[3]

                U = Un

                cinetique_all[m][0].append(U[0])
                cinetique_all[m][1].append(U[1])

    return cinetique_all


def animation_traj_2D(liste_p, liste_cin ,v): 
    '''
    Liste des corps, liste cinétique (des positions) , vitesse d'affichage
 
    Permet de visualiser une animation des trajectoires après avoir calculé toutes les positions au préalable (PAS en faisant pas par pas) (EN 2D)
    '''
  
    lines = [plt.plot([], [], '.-', color = p.color)[0] for p in liste_p]
    plt.plot(soleil.x, soleil.y, 'o', color = soleil.color) 
    tps = len(liste_cin[0][0])
    nb_planetes = len(liste_p)

    maxi = 1.1*max(liste_cin[-1][0])

    fig = plt.figure(1) # initialise la figure
    
    plt.xlim(-maxi, maxi)
    plt.ylim(-maxi, maxi)

    def init():

        for line in lines:
            line.set_data([], [])

        return lines

    def animate(i): 

        for k in range(nb_planetes):
            
            line = lines[k]
            x = liste_cin[k][0][v*i]
            y = liste_cin[k][1][v*i]


            line.set_data(x, y)

        return lines
    
    ani = animation.FuncAnimation(fig, animate, frames=range(int(tps/v)), init_func=init, blit=True, interval=1, repeat=True)

    plt.show()



def animation_traj_Lune_2D(liste_p, liste_cin ,v): 
    '''
    Liste des corps, liste cinétique (des positions) , vitesse d'affichage
 
    Pareil que animation_traj_2D sauf qu'il y a ici un changement d'échelle pour mieux visualiser 
    '''
    lines = [plt.plot([], [], '.-', color = p.color)[0] for p in liste_p]
    plt.plot(soleil.x, soleil.y, 'o', color = soleil.color) 
    tps = len(liste_cin[0][0])
    nb_planetes = len(liste_p)

    maxi = 1.1*max(liste_cin[-1][0])
    Terre=liste_planetes[2]
    fig = plt.figure(1) # initialise la figure


    plt.xlim(1.3E11, 1.6E11)
    plt.ylim(0E10,3E10)

    def init():

        for line in lines:
            line.set_data([], [])

        return lines

    def animate(i): 

        for k in range(nb_planetes):
            
            line = lines[k]
            x = liste_cin[k][0][v*i]
            y = liste_cin[k][1][v*i]


            line.set_data(x, y)

        return lines
    
    ani = animation.FuncAnimation(fig, animate, frames=range(int(tps/v)), init_func=init, blit=True, interval=1, repeat=True)

    plt.show()


# ANIMATION DE LA LUNE 2D
Lune = corps(2,0.00257,0.0549,7.36E22,'Lune',"y",systeme_solaire)
Lune.xi = liste_planetes[2].perihelie+356410E3
Lune.yi = 0
Lune.vxi = 0
Lune.vyi = 1.03756*liste_planetes[2].vyi
Lune.x =liste_planetes[2].perihelie+356410E3
Lune.y = 0
Lune.vx = 0
Lune.vy = 1.03756*liste_planetes[2].vyi

def simulation_terre_lune_soleil(beta, temps,dt):
    '''
    Fonction permettant de retourner les listes de trajectoires du système Terre-Lune-Soleil (ensuite on exécute animation_traj_2D(liste_p, liste_cin ,v) si l'on veut visualier l'animation
    '''
    liste=[Lune, liste_planetes[2]]
    for p in liste:
        p.refresh_cin()
        
    init_system(liste)
    Lune.xi = (liste_planetes[2].perihelie+356410E3)
    Lune.yi = 0
    Lune.vxi = 0
    Lune.vyi = 1.03756*liste_planetes[2].vyi
    Lune.x =(liste_planetes[2].perihelie+356410E3)
    Lune.y = 0
    Lune.vx = 0
    Lune.vy = 1.03756*liste_planetes[2].vyi
    return liste, systeme_precalc_fast(liste, beta, temps, dt)

beta = 2
liste_lune, liste_cin_lune = simulation_terre_lune_soleil(beta, 1E6, 1E2)


## ANIMATION EN VERSION 3D

def systeme_temps_reel_3D(liste_corps, beta, tmax, dt, freq):
    '''
    Calcul et affiche,grâce à Euler, pas par pas, la trajectoire des planètes en 3D
    '''
    
    nb_planetes = len(liste_corps)
    fref = force(2, astre_ref, astre_ref, systeme_solaire)
    k = fref.k
    temps = 0
    max = liste_corps[-1].a
    fig = plt.figure()
    plt.ion()

    c = 1
    while temps < tmax:
        
        ax = fig.add_subplot(projection='3d')
    
        for p in range(nb_planetes):

            planete = liste_corps[p]
            fs = force(beta, soleil, planete, systeme_solaire)
            planete.forces[0] = fs

            for k in range(p+1, nb_planetes, 1):

                planete2 = liste_corps[k]
                f12 = force(beta, planete, planete2, systeme_solaire)
                f21 = force(beta, planete2, planete, systeme_solaire)

                planete.forces[k] = f21
                planete2.forces[p+1] = f12

        for p in liste_corps:
            p.resultante = force_resultante(p)

        for p in liste_corps:
            p.vx = p.vx + (p.resultante[0]/p.mass)*dt
            p.vy = p.vy + (p.resultante[1]/p.mass)*dt

            p.x = p.x + p.vx*dt
            p.y = p.y + p.vy*dt
        
            ax.scatter(p.x,p.y,0,'.',color = p.color)
        ax.scatter(soleil.x, soleil.y, 0, 'o', color = soleil.color)

       
       #POUR AFFICHER L'EVOLUTION DU NOMBRE DE JOURS
       
        ax.axis([-2E11,2E11,-1E11,1E11])

        nbre=temps//(3600*24)
        if nbre>=c and nbre <c+1:
            #print(nbre, liste_planetes[2].vecteur_directeur(Lune))
            c += 1
        ann=ax.text2D(0, 0,'%d jours' % (nbre),transform=ax.transAxes)


        plt.pause(1/freq)
        plt.show()
        plt.clf()
        temps+=dt
        

def simulation_temps_reel_3D(n, beta, temps, dt, freq):
    '''
    Pour un groupe de planète (soit les 4 planètes intérieures, soit les 4 planètes extérieure), exécute systeme_temps_reel_3D(liste, beta, temps, dt, freq)
    '''

    if n == 0:
        liste = liste_planetes[0:4]
    else:
        liste = liste_planetes[4:9]
        
    for p in liste:
        p.refresh_cin()
        
    init_system(liste)
    systeme_temps_reel_3D(liste, beta, temps, dt, freq)


#Définitions des données nécessaires pour la simulation du système contenant la Lune
Lune = corps(2,0.00257,0.0549,7.36E22,'Lune',"y",systeme_solaire)
Lune.xi = liste_planetes[2].perihelie+356410E3
Lune.yi = 0
Lune.vxi = 0
Lune.vyi = 1.03756*liste_planetes[2].vyi
Lune.x =liste_planetes[2].perihelie+356410E3
Lune.y = 0
Lune.vx = 0
Lune.vy = 1.03756*liste_planetes[2].vyi


def simulation_terre_lune_soleil_3D(beta, temps,dt,freq):
    
    '''Permet de visualiser le système Terre-Lune-Soleil grâce à une fonction systeme_temps_reel_3D'''
    liste=[Lune, liste_planetes[2]]
    for p in liste:
        p.refresh_cin()
        
    init_system(liste)
    Lune.xi = (liste_planetes[2].perihelie+356410E3)
    Lune.yi = 0
    Lune.vxi = 0
    Lune.vyi = 1.03756*liste_planetes[2].vyi
    Lune.x =(liste_planetes[2].perihelie+356410E3)
    Lune.y = 0
    Lune.vx = 0
    Lune.vy = 1.03756*liste_planetes[2].vyi
    systeme_temps_reel_3D(liste, beta, temps, dt, freq)



### KEPLER

def calcul_periode(X, Y, dX, dY, temps):
    '''Calcule période en fonction de X, Y et des pas dx et dy, renvoie une période de révolution si il y a une révolution, sinon la plage de temps '''
    #Quand vx passe de positive à négative
    i = 1
    while (not (dX[i] > 0 and dX[i+1] < 0)) and (i<len(temps)):
        i += 1

    return temps[i] - temps[0]


def Kepler():
    '''Permet de verifier la troisième loi de Kepler
    '''
    #création des listes nécesssaires

    list_period_theo=[]
    list_period_exp=[]
    list_loiKepler_exp=[]
    list_loiKepler_theo=[]
    ecart=[]

    nb_planetes = len(liste_planetes_circ)
    for k in range(nb_planetes):

        planet, temps, X, Y, dX, dY = revolution(k, 2, 5, 1E5)
        period_exp = calcul_periode(X, Y, dX, dY, temps)/(3600*24)#calculs des périodes expérimentales
        list_period_exp.append(period_exp)

        period_th = liste_planetes_circ[k].period/(3600*24)#calculs des périodes théoriques
        list_period_theo.append(period_th)


    print("Liste des valeurs théorique des périodes:\n",list_period_theo, "\n\n""Liste des valeurs expérimentales des périodes:\n",list_period_exp )

    for i in range(nb_planetes):        
        #On regarde si on a une erreur plus grande que 1 jour entre les valeurs expérimentales et théoriques

        if abs(list_period_theo[i]-list_period_exp[i])>1:
            print("\nOn a un écart de plus de 1 jour entre les périodes de révolution trouvées expérimentalement et théoriquement")

    print("\nOn trouve les mêmes périodes de révolution expérimpentalement et théoriquement pour toutes les planètes, à 1 jour près\n")

    for i in range(nb_planetes): #calculs des rapports T**2/a**3

        list_loiKepler_exp.append((list_period_exp[i])**2/(liste_planetes[i].a)**3)


    for i in range(nb_planetes):

        list_loiKepler_theo.append((list_period_theo[i])**2/(liste_planetes[i].a)**3)

    for ind in range(nb_planetes):
        ecart.append(list_loiKepler_exp[ind]/list_loiKepler_theo[ind])

    #Trace la courbe

    fig, ax = plt.subplots(figsize=(5.5,6))
    ax.plot([p.name for p in liste_planetes_circ], ecart, 'x', color='g')
    ax.plot([p.name for p in liste_planetes_circ], [1 for ec in ecart], color='r')

    plt.title('Validité de la Loi de Kepler')
    ax.set_ylabel(ylabel='Rapport entre Loi de Kepler expérimentale et théorique', color='b', size=10)
    plt.ylim([0.7,1.25])
    plt.yticks(np.arange(0.7, 1.3, 0.05))

    for tick in ax.get_xticklabels():
        tick.set_rotation(45)
    plt.show()
    
    
### ETUDE ENERGIESs

def energie(planet, beta, X, Y, dX, dY, temps):
    '''Renvoie rayons, liste de temps, liste Ec, Ep et Em, et écarts entre énergies th et exp entre planète et soleil'''
    a = planet.a
    mass = planet.mass

    Eth = -k*mass/(2*a)

    E=[]
    #Liste des énergies mécaniques calculées à chaque temps et à chaque rayon pour une planète donnée
    R=[]
    #la liste des distances entre la planète considérée et le Soleil
    EC = []
    EP = []
    Ecarts = []
    #Listes des écarts entre les valeurs théoriques et expérimentales de l'énergie totale d'une planète

    nb_points = len(X)

    for i in range(nb_points):

        r = m.sqrt(X[i]**2+Y[i]**2)
        Ep=k*mass/((1-beta)*r**(beta-1))
        Ec=k*mass/(2*r**(beta-1))
        Etot=Ep+Ec

        E.append(Etot)
        R.append(r)
        EC.append(Ec)
        EP.append(Ep)
        Ecarts.append(abs(Etot - Eth))

    return (E, R, EC, EP, temps, Ecarts, Eth)

def energie_meca(planet, beta, X, Y, temps, tau):
    '''
    Affiche l'énergie mécanique en fonction du rayon(pour voir si c'est bien elliptique) et en fonction du temps (pour voir s'il y a pas divergence)
    '''
    a = planet.a
    mass = planet.mass

    Eth = -k*mass/(2*a)

    E=[]
    R=[]
    Ecarts = []

    nb_points = len(X)

    for i in range(nb_points):

        r = m.sqrt(X[i]**2+Y[i]**2)
        Ep=k*mass/((1-beta)*r**(beta-1))
        Ec=k*mass/(2*r**(beta-1))
        Etot=Ep+Ec

        E.append(Etot)
        R.append(r)
        Ecarts.append(abs(Etot - Eth))

    fig, (ax1, ax2) = plt.subplots(2, sharex=False, sharey=True)
    fig.suptitle("Energie mécanique tau=" + str(tau))
    ax1.plot(R, E, '.', R, [Eth for r in R], '-')
    ax1.set_ylabel('Energie mécanique')
    ax2.plot(temps, E, '.', temps, [Eth for t in temps], '-')
    ax2.set_ylabel('Energie mécanique')
    ax2.set_xlabel('Temps')
    ax2.set_title('Rayon')


    for ax in fig.get_axes():
        ax.label_outer()
    plt.show()

def etude_energie(planet, E, R, EC, EP, Ecarts, Eth, temps):
    '''Version graphique de la fonction énergie'''

    a = planet.a
    b = planet.b


    figb, (ax1b, ax2b, ax3b) = plt.subplots(3, sharex=True, sharey=True)
    figb.suptitle("Energie en fonction du temps",fontsize=9)
    ax1b.plot(temps, EC,'-' )
    ax1b.set_title("Energie cinétique",fontsize=7)
    ax2b.plot(temps, EP, '-')
    ax2b.set_title("Energie potentielle",fontsize=7)
    ax3b.plot(temps, E, '-')
    ax3b.set_title("Energie mécanique",fontsize=7)

    for ax in figb.get_axes():
        ax.label_outer()
    plt.show()

    fig, (ax1, ax2, ax3) = plt.subplots(3, sharex=True, sharey=True)
    fig.suptitle("Energie en fonction de la distance à l'étoile",fontsize=9)
    ax1.plot(R, EC,'-' )
    ax1.set_title("Energie Cinétique",fontsize=7)
    ax2.plot(R, EP, '-')
    ax2.set_title("Energie Potentielle",fontsize=7)
    ax3.plot(R, E, '-')
    ax3.set_title("Energie Mécanique",fontsize=7)
    
    for ax in fig.get_axes():
        ax.label_outer()
    plt.show()
    
    fig2, (ax1, ax2) = plt.subplots(2, sharex=True, sharey=False)
    ax1.set_title("Energie en fonction de la distance à l'étoile",fontsize=9)
    ax1.plot(R,E,'-', R, EC, '-', R, EP, '-') 
    ax1.plot(R, [Eth for r in R], '-')
    ax2.set_title("Ecart à la valeur théorique d'énergie")
    ax2.plot(R, Ecarts)
    ax2.plot(R, [abs(Eth) for r in R], '-')
    
        
    for ax in fig.get_axes():
        ax.label_outer()
    plt.show()


### ETUDE DE TAU / PAS D'INTEGRATION

def etude_tau_odeint(id_planet, beta, nb_periodes):
    '''
    Pour une planète, calcule la trajectoire avec odeint. Fait plusieurs tracé en variant le pas d'intégration.
    '''
    DTAU = [pow(10,i) for i in range(4,8,1)]
    K = [revolution(id_planet, beta, nb_periodes, dt) for dt in DTAU]
    for i in range(len(DTAU)):
        
        dt = DTAU[i]
        k = K[i]
        
        print_trajectoire_dtauodeint(k[0], beta, k[2], k[3], "Trajectoire de la Terre / dt = " + str(dt))


def etude_tau_Euler(id_planet, beta, tmax):
    '''
    Pour une planète, calcule avec l'intégration pas par pas avec Euler. On trace l'énergie mécanique pour différent pas d'intégration
    Donc permet de montrer l'influence du pas d'intégration sur l'énergie (AVEC EULER)

    '''
    DTAU = [pow(10,i) for i in range(2,7,1)]
    planet = liste_planetes[id_planet]
    liste = [planet]
    K = []
    for dtau in DTAU:
        temps = np.linspace(0, tmax, int(tmax/dtau))
        init_system(liste)
        [[X, Y]] = systeme_precalc_fast(liste, beta, tmax, dtau)
        energie_meca(planet, beta, X, Y, temps, dtau)
        



def etude_tau_RK4(id_planet, beta, tmax):
    '''
    Pour une planète, calcule avec l'intégration pas par pas avec RK4. On trace l'énergie mécanique pour différent pas d'intégration
    Donc permet de montrer l'influence du pas d'intégration sur l'énergie (AVEC RK4)
    '''
    DTAU = [pow(10,i) for i in range(2,7,1)]
    planet = liste_planetes[id_planet]
    liste = [planet]
    K = []
    for dtau in DTAU:
        temps = np.linspace(0, tmax, int(tmax/dtau))
        init_system(liste)
        [[X, Y]] = systeme_precalc_RK4(liste, beta, tmax, dtau)
        energie_meca(planet, beta, X, Y, temps, dtau)
        
        
        
### COMMANDES UTILES A EXECUTER

beta = 2

## ETUDE D'UNE PLANETE

simulation(2, 2, 5, 1E4)   #Pour afficher la trajectoire d'une planète (ici 2 donc pour la Terre) et les tracés énergétique 


#Kepler()   #Pour afficher les résultats montrant la validation de la troisième loi de Kepler


## ETUDE DE PLUSIEURS PLANETES QUI INTERAGISSENT ENTRE ELLES

#AFFICHAGE APRES TOUS LES CALCULS

# liste = liste_planetes[:4]        #Ici on regarde l'étude pour les 4 premières planètes
# init_system(liste)                #Pour initialiser les forces subies pour chaque planètes
# test = systeme_precalc_RK4(liste, beta, 1E7, 1E3)    #Calcul les trajectoires pour toutes les planètes considérées avec la méthode de RK4
# animation_traj_2D(liste, test, 10)                   #Montre l'animation découlant des trajectoires calculées précédemment
# print_trajectoires(liste, beta, test)



#OU AVEC AFFICHAGE ET CALCUL EN MEME TEMPS

#Affiche en 2D la trajectoire des 4 premières planètes
# simulation_temps_reel_2D(0, 2, 1E6, 1E4, 10000)


#Affiche en 3D la trajectoire des 4 premières planètes

#simulation_temps_reel_3D(0, 2, 1E6, 1E4, 10000)



##POUR OBSERVER LA DIVERGENCE SELON LA VALEUR DE DTAU POUR MERCURE ET URANUS

# init_system([liste_planetes[1],liste_planetes[6]]) 
# liste_cin = systeme_precalc_fast([liste_planetes[1],liste_planetes[6]], 2, 3E8, 1E5, tmin=0)
# print_trajectoires([liste_planetes[1],liste_planetes[6]], 2, liste_cin)
#  
 
## LUNE
#Permet d'afficher la trajectoire de la Lune et de la Terre autour du Soleil

#NOTES: On peut visualiser l'animation de la Lune:
#Soit en 3D pas par pas
#Soit en 2D avec animation après avoir tous calculé


#EN 2D
# liste_lune, liste_cin_lune = simulation_terre_lune_soleil(beta, 1E6, 1E2)
# animation_traj_Lune_2D(liste_lune, liste_cin_lune, 10) 

#EN 3D
# simulation_terre_lune_soleil_3D(2, 1E8, 1E5, 10000) 

##ETUDE DE TAU

# etude_tau_Euler(2, 2, 5E7)       #pour visualiser les courbes montrant l'energie mécanique en fonction du rayon et du temps pour la planète d'indice 2 (donc la Terre. Permet de montrer l'influence du pas d'intégration sur l'énergie  (AVEC EULER)

# etude_tau_odeint(2, 2, 5)     #Affiche la trajectoire de la planète 2 (donc la Terre) calculée pour différent dtau avec Odeint

# etude_tau_RK4(2, 2, 5E7)  #pour visualiser les courbes montrant l'energie mécanique en fonction du rayon et du temps pour la planète d'indice 2 (donc la Terre)