# -*- coding: utf-8 -*-
"""
Modèle décrivant la répartition des nématodes C.Elegans en deux dimension.
"""

# =============================================================================
# Bibliothèques
# =============================================================================
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# =============================================================================
# Fonctions                               ( opérateur * pour AND, + pour OR )        
# =============================================================================


# Fonction d'initialisation  : 
def initialisation(nb_sites_horizontal, nb_sites_vertical, nb_alleles) :
    nb_sites = nb_sites_vertical*nb_sites_horizontal
    # Initialisation de la matrice (que des zéros)
    matrice0 = np.zeros((nb_sites_vertical, nb_sites_horizontal, nb_alleles+1))  
    # Positions des sources
    positions_sources = np.random.binomial(1, proba_nourriture, nb_sites)  
    # Score des sources (entre -29 et 3)
    ages_sources = positions_sources*np.random.randint(score_pas_de_source+1, score_nouvelle_source+1, nb_sites)
    # Mise à -30 des autres sources
    ages_sources[positions_sources==0] = score_pas_de_source
    # Intégration dans le dernier élément des uplets de la matrice 
    matrice0[:,:,nb_alleles] = np.reshape(ages_sources, (nb_sites_vertical, nb_sites_horizontal))
    
    # Position des allèles (moit'-moit') sur les sources
    moitie_1 = range(0, int(nb_sites_horizontal/2))
    moitie_2 = range(int(nb_sites_horizontal/2), nb_sites_horizontal)
    matrice0[:, moitie_1, 0] = np.reshape(positions_sources, (nb_sites_vertical, nb_sites_horizontal))[:, moitie_1]
    matrice0[:, moitie_2, 1] = np.reshape(positions_sources, (nb_sites_vertical, nb_sites_horizontal))[:, moitie_2]
    return(matrice0)



   
# Fonction d'apparition des nouvelles sources
def nouvelles_sources(matrice):
    # Positions des nouvelles sources 
    nouvelles_sources = np.random.binomial(1, proba_nourriture, nb_sites_vertical*nb_sites_horizontal)
    nouvelles_sources = np.reshape(nouvelles_sources, (nb_sites_vertical, nb_sites_horizontal))
    # Emplacements sans sources (de score -30)
    emplacements_vides = (matrice[:, :, nb_alleles] == score_pas_de_source)
    
    # Emplacements des nouvelles sources dans les places vides, fois 4 (score nouvelle source non colonisée = score nouvelle source + 1)
    matrice[:,:,nb_alleles] = matrice[:,:,nb_alleles] + (nouvelles_sources*emplacements_vides)*(score_nouvelle_source-score_pas_de_source+1)
    return(matrice)
 


# Fonction de colonisation des sources
def colonisation(matrice):
    
   # sources non-occupées par des nématodes (sites dont les fréquences = 0 et dont le score = score nouvelle source + 1)
   sources_non_decouvertes = (matrice[:,:,nb_alleles]==score_nouvelle_source+1) 
   # sources en migration (entre -1 et -29)
   sources_en_migration = (matrice[:,:,nb_alleles]<score_migration)*(matrice[:,:,nb_alleles]>score_pas_de_source)

   # S'il y a des sources en migration ET des sources non découvertes...
   if np.any(sources_en_migration)*np.any(sources_non_decouvertes) :
       # En énumérant les sources non découvertes
       for k in np.transpose(np.where(sources_non_decouvertes)) :
           k[0]=int(k[0])
           k[1]=int(k[1])
           sources_colonisatrices = []
           
           # En énumérant les sources en migration se situant dans le rayon de migration
           borne_min_vertical = max(k[0]+score_pas_de_source, k[0]-rayon_migration, 0)  
           borne_max_vertical = min(k[0]-score_pas_de_source, k[0]+rayon_migration, nb_sites_vertical-1)  
           borne_min_horizontal = max(k[1]+score_pas_de_source, k[1]-rayon_migration, 0)
           borne_max_horizontal = min(k[1]-score_pas_de_source, k[1]+rayon_migration, nb_sites_horizontal-1)
           
           
           for j in np.transpose(np.where(sources_en_migration[borne_min_vertical:borne_max_vertical+1,borne_min_horizontal:borne_max_horizontal+1])) :  # attention indice fin non compris
               j[0] = borne_min_vertical + int(j[0])
               j[1] = borne_min_horizontal + int(j[1])
               
               # si la distance <= score migration et si elle ne dépasse pas le rayon de migration
               if abs(j[0]-k[0]) <= -matrice[j[0], j[1], nb_alleles] and abs(j[1]-k[1]) <= -matrice[j[0], j[1], nb_alleles] :
                   sources_colonisatrices = np.append(sources_colonisatrices, j)
                   sources_colonisatrices = sources_colonisatrices.astype(int)                     
                       
           sources_colonisatrices = np.reshape(sources_colonisatrices, (int(len(sources_colonisatrices)/2),2))  
           
           if len(sources_colonisatrices) != 0 :
               # Fréquences sur l'ensemble des sources colonisatrices (si pas de sources colonisatrices, ce vecteur = (0,0))        
               vecteur_allele = np.sum(matrice[sources_colonisatrices[:,0], sources_colonisatrices[:,1], 0:nb_alleles], axis=0)
             
               # Méthode déterministe
               if type_colonisation == 'deterministe' : 
                   matrice[k[0],k[1],0:nb_alleles] = vecteur_allele/np.sum(vecteur_allele)
                   
               # Avec une loi de poisson
               if type_colonisation == 'stochastique' : 
                   # Nombre de parents de la nouvelle source (on en veut au moins un !)
                   nb_parents = np.random.poisson(lamb) 
                   while nb_parents == 0 : nb_parents = np.random.poisson(lamb)         
                   # Choix des allèles de ces parents dans la jauge totale de toutes les sources colonisatrices (loi multinomiale)    
                   # Le vecteur alleles_parents contient à chaque indice, le nombre de parents portant l'allèle correspondant à l'indice.
                   alleles_parents = np.random.multinomial(nb_parents, vecteur_allele/np.sum(vecteur_allele))    
                   # Les proportions dans la sous-population des parents devient la jauge de la nouvelle source
                   matrice[k[0],k[1],0:nb_alleles] = alleles_parents/np.sum(alleles_parents)
                   # Le score de la source colonisée devient 3.5 
               matrice[k[0],k[1],nb_alleles] = score_nouvelle_source+0.5           

   return(matrice)




# Fonction temps
def pas_de_temps(matrice):
   
   if jauges_partout == 'non' :
        # On enlève les jauges des sources dont le score est de -19 (qui vont passer à -20 et disparaître)
        sources_disparues = (matrice[:,:,nb_alleles]==score_pas_de_source+1)
        if np.sum(sources_disparues) != 0 :
            matrice[sources_disparues,0:nb_alleles]=np.zeros((np.sum(sources_disparues),2))
  
   # sources de score 4 occupées par des nématodes
   sources_decouvertes = (matrice[:,:,nb_alleles]==score_nouvelle_source+0.5)   
   # sources ayant des scores entre 3 et -29 
   sources_actives = (matrice[:,:,nb_alleles]>score_pas_de_source)*(matrice[:,:,nb_alleles]<=score_nouvelle_source)
    
   # Diminution des scores de 1
   matrice[sources_decouvertes, nb_alleles]=matrice[sources_decouvertes, nb_alleles]-0.5 
   matrice[sources_actives, nb_alleles]=matrice[sources_actives, nb_alleles]-1 
   
   # Gestion des sources non colonisées : certaines restent, d'autres disparaissent (avec une proba q à chaque pas de temps)
   sources_non_colonisees = (matrice[:,:,nb_alleles]==score_nouvelle_source+1)
   # Ajout de -34 au score des sources qui vont disparaître (4 -> -30). Ces sources sont tirées selon une proba q (proba_disparition_nourriture).
   matrice[sources_non_colonisees, nb_alleles] = matrice[sources_non_colonisees,nb_alleles] + np.random.binomial(1,proba_disparition_nourriture, np.sum(sources_non_colonisees))*(score_pas_de_source-score_nouvelle_source-1) 
 
   return(matrice)
   
      
   
   
# Fontion d'évolution en fonction du temps (très détaillée)
def evolution_detail(nb_sites_horizontal, nb_sites_vertical, nb_alleles, temps):
    matrice = initialisation(nb_sites_horizontal, nb_sites_vertical, nb_alleles)
    print ("Environnement au temps initial")
    print(matrice)
    couleur(matrice, nb_alleles)
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
            couleur(matrice, nb_alleles)
    if affiche_tous_les_histogrammes=='non' :
        couleur(matrice, nb_alleles)
    print ("Environnement au temps final")
    return(np.around(matrice, decimals=2))
    
    
    
    
# Fontion d'évolution en fonction du temps (juste histogrammes)
def evolution(nb_sites_horizontal, nb_sites_vertical, nb_alleles, temps):
    matrice = initialisation(nb_sites_horizontal, nb_sites_vertical, nb_alleles)
    couleur(matrice, nb_alleles)
    print ("Environnement au temps initial")
    for i in range(0,temps) :
        matrice = nouvelles_sources(matrice)
        matrice = colonisation(matrice)   
        matrice = pas_de_temps(matrice)
        if affiche_tous_les_histogrammes=='oui' :
            couleur(matrice, nb_alleles)
            print ("Environnement au temps %i" %i)
    if affiche_tous_les_histogrammes=='non':
        couleur(matrice, nb_alleles)
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
    nb_sites = nb_sites_vertical*nb_sites_horizontal
    # Position des nouvelles sources pour la nouvelle année
    apparition_sources = np.reshape(np.random.binomial(1, proba_nourriture, nb_sites), (nb_sites_vertical, nb_sites_horizontal))  # Donne les emplacements où nouvelle nourriture
    jauges_existantes = (matrice[:,:,nb_alleles]<=score_nouvelle_source)*(matrice[:,:,nb_alleles]>score_pas_de_source)            # Donne les emplacement où une jauge était présente
    positions_nouvelles_sources = apparition_sources*jauges_existantes
    # Age des nouvelles sources mis à jour directement sur la matrice
    matrice[:,:,nb_alleles] = positions_nouvelles_sources*np.reshape(np.random.randint(score_migration, score_nouvelle_source+1, nb_sites), (nb_sites_vertical, nb_sites_horizontal))
    # On efface la jauge sur les sites sans nouvelle source
    matrice[positions_nouvelles_sources==0,0:nb_alleles] = 0
    # Mise à -30 des autres sources
    matrice[positions_nouvelles_sources==0, nb_alleles] = score_pas_de_source 
    
    if jauge == 'unie' : 
        # On tire au hasard l'allèle qui occupe les nouvelles sources en fonction de la jauge précédente
        matrice[positions_nouvelles_sources==1,0] = np.random.binomial(1,matrice[positions_nouvelles_sources==1,0],len(matrice[positions_nouvelles_sources==1,0]))
        matrice[positions_nouvelles_sources==1,1] = 1-matrice[positions_nouvelles_sources==1,0]
    return(matrice)



# Fonction d'évolution en fonction du temps et de la matrice initiale fournie à chaque printemps (ne trace rien)
def evolution_itere(matrice):
    for i in range(0,temps) :
        matrice = nouvelles_sources(matrice)
        matrice = colonisation(matrice)
        matrice = pas_de_temps(matrice)
    return(np.around(matrice, decimals=4))



# Dessine l'histogramme avant et après chaque hiver
def evolution_annees(nb_sites_horizontal, nb_sites_vertical, nb_alleles, nb_annees, temps_annees, jauge):
    matrice = initialisation(nb_sites_horizontal, nb_sites_vertical,  nb_alleles)
    couleur(matrice, nb_alleles)
    print("Environnement au début de l'année 0")
    matrice = evolution_itere(matrice)
    couleur(matrice, nb_alleles)
    print ("Environnement à la fin de l'année 0")
    for i in range(1,nb_annees+1):
        matrice = matrice_printemps(matrice, jauge)
        couleur(matrice, nb_alleles)
        print("Environnement au début de l'année %i" %i)
        matrice = evolution_itere(matrice)
        couleur(matrice, nb_alleles)
        print ("Environnement à la fin de l'année %i" %i)
    return(matrice)



    


# =============================================================================
# Illustration
# =============================================================================
    
    
# Fontion pour tracer les matrices colorées
def couleur(matrice, nb_alleles):
    matrice_graph = matrice
    matrice_graph[matrice[:,:,0]+matrice_graph[:,:,1]==0,0:nb_alleles]=-1
    plt.figure()
    for i in range(0,nb_alleles) :
        plt.subplot(1,nb_alleles,i+1)
        sns.heatmap(matrice_graph[:,:,i], cmap="YlGnBu", cbar=False)          # cmap="PiYG"
        plt.title('Allèle %i' %(i+1))
    plt.show()


# =============================================================================
# Etude de l'évolution de la répartition des allèles en temps
# =============================================================================

nb_sites_vertical = 50       # Nombres de sites en longueur (verger)            /!\ programmation : vertical = ligne
nb_sites_horizontal = 50     # Nombres de sites en largeur (verger)                                 horizontal = colonne
nb_alleles = 2               # Nombres d'allèles étudiés

proba_nourriture = 0.2               # Proba d'pparition de la nourriture   
proba_disparition_nourriture = 0.05  # Proba qu'une source non colonisée disparaisse   

score_nouvelle_source = 3     # Score des sources 
score_migration = 0
score_pas_de_source = -30

rayon_migration = 5       # Distance maximale qu'un vers peut atteindre en migrant à partir de sa colonie
 
temps = 100               # Intervalle de temps étudié

type_colonisation = 'stochastique'   # Methode pour coloniser une nouvelle jauge : 'deterministe' = moyenne des fréquence, 'stochastique' = choix des parents par loi de poisson
lamb = 4                            # Dans le cas : type_colonisation = 'stochastique', paramètre de la loi de poisson qui choisit le nb de parents


jauges_partout = 'non'                     # Afficher des jauges partout (garde en mémoire la dernière population ayant occupé le site)
affiche_tous_les_histogrammes = 'non'      # Afficher les histogrammes à chaque pas de temps.  


nb_annees  = 5
temps_annees = 120
jauge = 'unie'             # 'unie' : A chaque printemps, on ne prends en compte qu'un individu par jauges encore présentes.
                           # 'complete' :  A chaque printemps, on ne prends en compte que les jauges encore présentes.







# Pour avoir les histogrammes uniquement.
matrice_1 = evolution(nb_sites_horizontal, nb_sites_vertical, nb_alleles, temps)   
# Pour voir l'évolution sur plusieurs années (entre-coupées d'hivers)
matrice_2 = evolution_annees(nb_sites_horizontal, nb_sites_vertical, nb_alleles, nb_annees, temps_annees, jauge)



# VERIFICATION
# Pour avoir toutes les matrices après chaque action.
# evolution_detail(nb_sites_horizontal, nb_sites_vertical, nb_alleles, temps)






