FORCAGE STOCHASTIQUE DE LA TURBULENCE
===

Ce développement implémente un forçage stochastique de la turbulence. Ce forçage est calculé tous les *dtfst* pas de temps (en seconde) et est interpolé entre ces deux calculs consécutifs. 


La fonction de forçage est déninie de la façon suivante:

$$
\begin{align}
F_{n\Delta t}(i,j,k) = F_{T}.W(i,j).F^{-1}[\sum_{n} N_{n\Delta t}(kx,ky,n).\psi(kx,ky,n).\phi_{n}(k)]
\end{align}
$$

$ F_{T} $ est une constante (qui se déclinera ensuite en $ F_{u}, F_{v} $ ) pour forcer la densité

$ W(i,j) $ est une fonction de forme du forçage (pour l'instant seulement en y) dans l'espace physique

$$ W(i,j) = e^{-((yr(i,j)-y_{F})/(\Delta y_{F}))^{2}} $$

$N_{n\Delta t}(kx,ky,n)$ est un bruit blanc. C'est une fonction aléatoire reproductible dans le temps.

$ \psi(kx,ky,n) $ définit l'amplitude du signal

$$ 
\psi (kx,ky,n)= \underbrace{$a_{n}}_{mode} \underbrace{e^{-(kx^{2} + ky^{2} - kF^{2})/(\Delta kF^{2}) }}_{horizontal}
$$

$\phi_{n}$ sont les modes verticaux

Implémentation dans le code
===

- cppdefs.h : ajout d'une clé *FSTURB* pour définir le forçage
- ana_fsturb.F : nouveau module qui regroupe les fonctions d'initialisation et de calcul du forçage dans le temps
- main.F : initialisation du forçage ( *call init_fsturb*)
- step.F : calcul d'un nouveau forçage tous les $\Delta t$ pas de temps
- step3d_t.F : ajout du forçage interpolé en temps et en z
 
Paramètres
===

- *Nmodefst* : Nombre de modes
- dtfst = 1.*86400
- *afst* : $a_{n}$ amplitude des différents modes
- *xfstmid* :
- *xfstwid* :
- *yfstmid* : $ y_{F} $
- *yxfstwid* : $ \Delta y_{F} $
- *FTfst* : $ F_{T} $
- *kFfst* : $ kF $
- *dkFfst* : $ \Delta kF $

Variables
===

Toutes les variables relatives au forçage sont déclarées dans le module *fsturb*, dans le fichier *ana_fsturb.F*

- *Ffst(i,j,k,t)* : $ F_{n\Delta t}(i,j,k) $
- *FTfst* : $ F_{T} $ défini en parameter dans le module
- *Wfst* : W(i,j)
- *Noisefst* : $ N_{n\Delta t}(kx,ky,n) $ 
- *afst* : $a_{n}$ amplitude des différents modes
- *psifst* : $\psi(kx,ky,n)$
- *phir* : modes calculés dans vmodes

