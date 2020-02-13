# Run_MOZ

Le répertoire Run_MOZ doit se situer dans CONFIGS


Installation
=============
git clone https://github.com/slgentil/Run_MOZ.git

Compilation
============
Aller dans le répertoire CONFIGS/Run_MOZ 
Ce répertoire ne contient que les sources modifiés de la configuration courante par rapport à CROCO de base   
Lancer jobcomp  
L'exécutable croco et un lien symbolique vers l'exécutable xios_server.exe sont créés dans le répertoire courant

Lancement sur Datarmor
======================
_chain_datarmor.py workdir nbchain elaptim resolution jobname restart_  
- workdir : répertoire qui sera créé sous DATAWORK ou SCRATCH selon la variable WORK du script python
- nbchain : nombre de chainages  
- elaptim : temps elapsed pour chaque chainage HH:MM:SS  
- resolution : dépend de la configuration (256,512,1024)
- jobname : nom générique des batchs    
- restart : 0 (initial) or 1 (restart)  

Le répertoire workdir est créé dans DATAWORK ou SCRATCH  
Ce répertoire contient des sous-répertoires t1,t2,... , un répertoire par chainage  
Le script python prépare le batch et vous indique comment le lancer.  

