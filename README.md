CROCO
=====

Le répertoire croco contient le code CROCO de base  

Le répertoire CONFIGS contient les différentes configurations:  
- Run_MEDDY : Meddy dans une configuration idéalisée
- Run_MOZ   : Tourbillon dans la configuration réaliste du canal du Mozambique
- Run_JETN  : Jet et marée interne dans une configuration idéalisée  

Le répertoire util contient des utilitaires python  
- clean_croco.py : pour supprimer tous les xios_client* et xios_server* (sauf le 0) à partir du niveau de répertoire courant et dans les deux niveaux inférieurs  
- kill_datarmor.py : pour tuer sur datarmor tous les chainages d'une simulation en cours  
- restart_datarmor.py : pour relancer sur datarmor une simulation de plusieurs chainages  


Installation
=============
git clone https://github.com/slgentil/croco.git

Compilation
============
Aller dans le répertoire Run_XXX  
Ce répertoire ne contient que les sources modifiés de la configuration courante par rapport à CROCO de base   
Lancer jobcomp  
L'exécutable croco et un lien symbolique vers l'exécutable xios_server.exe sont créés dans le répertoire courant

Lancement sur Datarmor
======================
_chain_datarmor.py workdir nbchain elaptim resolution jobname restart_  
- workdir : répertoire qui sera créé sous DATAWORK ou SCRATCH selon la variable WORK du script python
- nbchain : nombre de chainages  
- elaptim : temps elapsed pour chaque chainage HH:MM:SS  
- resolution : dépend de la configuration (4,2,1 pour JETN / 128,256,512,1024 pour MEDDY et MOZ)
- jobname : nom générique des batchs    
- restart : 0 (initial) or 1 (restart)  

Le répertoire workdir est créé dans DATAWORK ou SCRATCH  
Ce répertoire contient des sous-répertoires t1,t2,... , un répertoire par chainage  
Le script python prépare le batch et vous indique comment le lancer.  

