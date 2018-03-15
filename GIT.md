
Reference:

https://git-scm.com/book/en/v2/Git-Branching-Remote-Branches

Very important command in order to know status (branch, modifications, ...):
```
git status
```

Bring data from remote server
```
git clone https://github.com/slgentil/croco.git 
git log --oneline --decorate --graph --all
```

How to commit changes to master:
```
git add modified_files
git commit -m 'Description of the modifications'
git push origin master 
```

How to commit to a new branch (mybranch):
```
git add modified_files
git checkout -b mybranch
git commit -m 'Description of the modifications'
git push origin mybranch
git branch -d mybranch
```

How to merge mybranch to master:
```
git checkout master
git pull origin master
git merge mybranch
git branch -d mybranch
```


Create local branch from remote ones
```
git checkout -b ap_changes origin/ap_changes
git checkout -b sf_changes origin/sf_changes
```

Now merge ap_changes into master
```
git merge ap_changes
(CONFLICT (content): Merge conflict in overview/plot_snapshot_2d_noproj.py)
vim overview/plot_snapshot_2d_noproj.py
git add overview/plot_snapshot_2d_noproj.py
git commit
```

Delete ap_changes branch
```
git branch -d ap_changes
```

Checkout a file from another branch
```
git checkout mybranch
git checkout otherbranch -- dev/file.py
```

