# -*- coding: utf-8 -*-
"""
Modèle décrivant la répartition des nématodes C.Elegans en une dimension.
"""

# =============================================================================
# Bibliothèques
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt


# =============================================================================
# Fonctions                               ( opérateur * pour AND, + pour OR )        
# =============================================================================


# Fonction d'initialisation  :
def initialisation(nb_sites, nb_alleles) :
    # Initialisation de la matrice (que des zéros)
    matrice0 = np.zeros((nb_alleles+1,nb_sites))  
    # Positions des sources
    positions_sources = np.random.binomial(1, proba_nourriture, nb_sites)  
    # Score des sources (entre -29 et 3)
    matrice0[nb_alleles] = positions_sources*np.random.randint(score_pas_de_source+1, score_nouvelle_source+1, nb_sites)
    # Mise à -30 des autres sources
    matrice0[nb_alleles, positions_sources==0] = score_pas_de_source
    
    # Position des allèles (moit'-moit') sur les sources
    moitie_1 = range(0, int(nb_sites/2))
    moitie_2 = range(int(nb_sites/2), nb_sites)
    matrice0[0, moitie_1] = positions_sources[moitie_1]
    matrice0[1, moitie_2] = positions_sources[moitie_2]
    return(matrice0)




   
# Fonction d'apparition des nouvelles sources
def nouvelles_sources(matrice):
    # Positions des nouvelles sources 
    nouvelles_sources = np.random.binomial(1, proba_nourriture, nb_sites)
    # Emplacements sans sources (de score -30)
    emplacements_vides = (matrice[nb_alleles,] == score_pas_de_source)
    
    # Emplacements des nouvelles sources dans les places vides, fois 4 (score nouvelle source non colonisée = score nouvelle source + 1)
    matrice[nb_alleles] = matrice[nb_alleles] + (nouvelles_sources*emplacements_vides)*(score_nouvelle_source-score_pas_de_source+1)
    return(matrice)
 
    
    
# Fonction de colonisation des sources
def colonisation(matrice):
    
   # sources non-occupées par des nématodes (sites dont les fréquences=0 et dont le score = score nouvelle source)
   sources_non_decouvertes = (matrice[nb_alleles,:]==score_nouvelle_source+1) 
   # sources en migration (entre -1 et -29)
   sources_en_migration = (matrice[nb_alleles,:]<score_migration)*(matrice[nb_alleles,:]>score_pas_de_source)
   
   # S'il y a des sources en migration ET des sources non découvertes...
   if np.any(sources_en_migration)*np.any(sources_non_decouvertes) :
       # En énumérant les sources non découvertes
       for k in np.transpose(np.where(sources_non_decouvertes)) :
           k = int(k)   
           sources_colonisatrices = []
           
           # En énumérant les sources en migration se situant dans le rayon de migration
           borne_min = max(k+score_pas_de_source, k-rayon_migration, 0)
           borne_max = min(k-score_pas_de_source, k+rayon_migration, nb_sites-1)  
           for j in np.transpose(np.where(sources_en_migration[borne_min:borne_max+1])) :  # attention indice fin non compris
               j = borne_min + int(j)
               
               # si la distance <= score migration et si elle ne dépasse pas le rayon de migration
               if abs(j-k) <= -matrice[nb_alleles,j] :
                   sources_colonisatrices = np.append(sources_colonisatrices, [j])
                   sources_colonisatrices = sources_colonisatrices.astype(int)
    
           # Fréquences sur l'ensemble des sources colonisatrices (si pas de sources colo, ce vect = (0,0))        
           vecteur_allele = np.sum(matrice[0:nb_alleles,sources_colonisatrices], axis=1)
           
           if sum(vecteur_allele) != 0 :
               # Rajout des fréquences sur la source colonisée
               
               # Méthode déterministe
               if type_colonisation == 'deterministe' : 
                   matrice[0:nb_alleles,k] = vecteur_allele/np.sum(vecteur_allele)
                   
               # Avec une loi de poisson
               if type_colonisation == 'stochastique' : 
                   # Nombre de parents de la nouvelle source (on en veut au moins un !)
                   nb_parents = np.random.poisson(lamb) 
                   while nb_parents == 0 : nb_parents = np.random.poisson(lamb)         
                   # Choix des allèles de ces parents dans la jauge totale de toutes les sources colonisatrices (loi multinomiale)    
                   # Le vecteur alleles_parents contient à chaque indice, le nombre de parents portant l'allèle correspondant à l'indice.
                   alleles_parents = np.random.multinomial(nb_parents, vecteur_allele/np.sum(vecteur_allele))    
                   # Les proportions dans la sous-population des parents devient la jauge de la nouvelle source
                   matrice[0:nb_alleles,k] = alleles_parents/np.sum(alleles_parents)
                   # Le score de la source colonisée devient 3.5 
               matrice[nb_alleles,k] = score_nouvelle_source+0.5           

   return(matrice)




# Fonction temps
def pas_de_temps(matrice):
   
   if jauges_partout == 'non' :
        # On enlève les jauges des sources dont le score est de -29 (qui vont passer à -30 et disparaître)
        sources_disparues = (matrice[nb_alleles,:]==score_pas_de_source+1)
        if sum(sources_disparues) != 0 :
            matrice[0:nb_alleles,sources_disparues]=np.zeros((nb_alleles,sum(sources_disparues)))
  
   # sources de score 4 occupées par des nématodes
   sources_decouvertes = (matrice[nb_alleles,:]==score_nouvelle_source+0.5)   
   # sources ayant des scores entre 3 et -29
   sources_actives = (matrice[nb_alleles,:]>score_pas_de_source)*(matrice[nb_alleles,:]<=score_nouvelle_source)
    
   # Diminution des scores de 1
   matrice[nb_alleles,sources_decouvertes]=matrice[nb_alleles,sources_decouvertes]-0.5 
   matrice[nb_alleles,sources_actives]=matrice[nb_alleles,sources_actives]-1 
   
   # Gestion des sources non colonisées : certaines restent, d'autres disparaissent (avec une proba q à chaque pas de temps)
   sources_non_colonisees = (matrice[nb_alleles,:]==score_nouvelle_source+1)
   # Ajout de -34 au score des sources qui vont disparaître (4 -> -30). Ces sources sont tirées selon une proba q (proba_disparition_nourriture).
   matrice[nb_alleles, sources_non_colonisees] = matrice[nb_alleles,sources_non_colonisees] +   np.random.binomial(1,proba_disparition_nourriture, sum(sources_non_colonisees))*(score_pas_de_source-score_nouvelle_source-1) 
 
   return(matrice)
   
      
   
   
# Fontion d'évolution en fonction du temps (très détaillée)
def evolution_detail(nb_sites, nb_alleles, temps):
    matrice = initialisation(nb_sites, nb_alleles)
    print ("Environnement au temps initial")
    print(matrice)
    histogramme(matrice)
    for i in range(0,temps) :
        print ("Environnement au temps", i)
        matrice = nouvelles_sources(matrice)
        print("Nouvelles sources")
        print(matrice)
        matrice = colonisation(matrice)   
        print("Colonisation")
        print(matrice)
        matrice = pas_de_temps(matrice)
        print("Pas de temps")
        print(matrice)
        if affiche_tous_les_histogrammes=='oui' :
            histogramme(matrice)
    if affiche_tous_les_histogrammes=='non' :
        histogramme(matrice)
    print ("Environnement au temps final")
    return(np.around(matrice, decimals=2))
    
    
    
    
# Fontion d'évolution en fonction du temps (juste histogrammes)
def evolution(nb_sites, nb_alleles, temps):
    matrice = initialisation(nb_sites, nb_alleles)
    histogramme(matrice)
    print ("Environnement au temps initial")
    for i in range(0,temps) :
        matrice = nouvelles_sources(matrice)
        matrice = colonisation(matrice)   
        matrice = pas_de_temps(matrice)
        if affiche_tous_les_histogrammes=='oui' :
            histogramme(matrice)
            print ("Environnement au temps %i" %i)
    if affiche_tous_les_histogrammes=='non':
        histogramme(matrice)
    print ("Environnement au temps final")
    return(np.around(matrice, decimals=1))
    







# =============================================================================
# Simulation des hivers (une partie de la population meurt)
# =============================================================================
    
    
    
# Retourne la matrice initiale au printemps,  
# qui contient uniquement les jauges encore présentes si jauge = 'complete', 
# qui contient uniquement un individu pour chaque jauge encore présente si jauge = 'unie'
# /!\ ne fonctionne qu'avec deux allèles
def matrice_printemps(matrice, jauge):
    #Position des nouvelles sources pour la nouvelle année
    apparition_sources = np.random.binomial(1, proba_nourriture, nb_sites)  # Donne les emplacements où nouvelle nourriture
    jauge_existante = sum(matrice[0:nb_alleles,:]) != 0                     # Donne les emplacement où une jauge était présente
    positions_nouvelles_sources = apparition_sources*jauge_existante
    #Age des nouvelles sources mis à jour directement sur la matrice
    matrice[nb_alleles] = positions_nouvelles_sources*np.random.randint(score_migration, score_nouvelle_source+1, nb_sites)
    #On efface la jauge sur les sites sans nouvelle source
    matrice[0:nb_alleles, positions_nouvelles_sources==0] = 0
    # Mise à -20 des autres sources
    matrice[nb_alleles, positions_nouvelles_sources==0] = score_pas_de_source
    
    if jauge == 'unie' : 
        #On tire au hasard l'allèle qui occupe les nouvelles sources en fonction de la jauge précédente
        matrice[0, positions_nouvelles_sources==1] = np.random.binomial(1, matrice[0, positions_nouvelles_sources==1], len(matrice[0,positions_nouvelles_sources==1]))
        matrice[1, positions_nouvelles_sources==1] = 1-matrice[0,positions_nouvelles_sources==1]
    return(matrice)




# Fonction d'évolution en fonction du temps et de la matrice initiale fournie à chaque printemps (ne trace rien)
def evolution_itere(matrice):
    for i in range(0,temps) :
        matrice = nouvelles_sources(matrice)
        matrice = colonisation(matrice)
        matrice = pas_de_temps(matrice)
    return(np.around(matrice, decimals=4))




# Dessine l'histogramme avant et après chaque hiver
def evolution_annees(nb_sites, nb_alleles, nb_annees, temps_annees):
    matrice = initialisation(nb_sites, nb_alleles)
    histogramme(matrice)
    print("Environnement au début de l'année 0")
    matrice = evolution_itere(matrice)
    histogramme(matrice)
    print ("Environnement à la fin de l'année 0")
    for i in range(1,nb_annees+1):
        matrice = matrice_printemps(matrice, jauge)
        histogramme(matrice)
        print("Environnement au début de l'année %i" %i)
        matrice = evolution_itere(matrice)
        histogramme(matrice)
        print ("Environnement à la fin de l'année %i" %i)
    return(matrice)





# =============================================================================
# Illustration
# =============================================================================
    
    
# Fonction pour tracer les graphiques
def graphique(matrice):
    plt.figure()
    positions = np.arange(nb_sites)
    p1=plt.plot(positions, matrice[0,:], marker='o')
    p2=plt.plot(positions, matrice[1,:], marker='v')
    plt.title("Répartition des allèles")  # Problemes avec accents (plot_directive) !
    plt.legend([p1, p2], ["Allèle 1", "Allèle 2"])
    plt.show()
    
    
# Fontion pour tracer les histogrammes
def histogramme(matrice):
    plt.figure()
    positions = np.arange(nb_sites)
    bins = [x - 0.5 for x in range(0, nb_sites+1)]
    plt.hist([positions, positions], bins = bins, weights = [matrice[0,:]*100, matrice[1,:]*100],
            edgecolor = 'black', histtype = 'barstacked', label = ['allèle 1', 'allèle 2'])
    plt.ylabel('Pourcentage')
    plt.xlabel('Sites')
    plt.title('Répartition des allèles (en %)')
    plt.legend()
    plt.show()
     
    
  


# =============================================================================
# Etude de l'évolution de la répartition des allèles en temps
# =============================================================================

nb_sites = 100              # Nombres de sites dans l'environnement
nb_alleles = 2              # Nombres d'allèles étudiés

proba_nourriture = 0.2               # Proba d'pparition de la nourriture   
proba_disparition_nourriture = 0.1  # Proba qu'une source non colonisée disparaisse   

score_nouvelle_source = 3     # Score des sources 
score_migration = 0
score_pas_de_source = -30

rayon_migration = 5       # Distance maximale qu'un vers peut atteindre en migrant à partir de sa colonie
 
temps = 600               # Intervalle de temps étudié

type_colonisation = 'stochastique'   # Methode pour coloniser une nouvelle jauge : 'deterministe' = moyenne des fréquence, 'stochastique' = choix des parents par loi de poisson
lamb = 4                             # Dans le cas : type_colonisation = 'stochastique', paramètre de la loi de poisson qui choisit le nb de parents


jauges_partout = 'non'                     # Afficher des jauges partout (garde en mémoire la dernière population ayant occupé le site)
affiche_tous_les_histogrammes = 'non'      # Afficher les histogrammes à chaque pas de temps.  


nb_annees  = 5
temps_annees = 120
jauge = 'unie'             # 'unie' : A chaque printemps, on ne prends en compte qu'un individu par jauges encore présentes.
                               # 'complete' :  A chaque printemps, on ne prends en compte que les jauges encore présentes.




# Temps continu (sans prendre en compte l'hiver)
matrice_1 = evolution(nb_sites, nb_alleles, temps)   

# Pour voir l'évolution sur plusieurs années (entre-coupées d'hivers)
matrice_2 = evolution_annees(nb_sites, nb_alleles, nb_annees, temps_annees)


# VERIFICATION
# Pour avoir toutes les matrices après chaque action.
# evolution_detail(nb_sites, nb_alleles, temps)    

