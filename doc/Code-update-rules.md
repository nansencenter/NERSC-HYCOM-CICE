# GIT-hub updating rules for the Ocean Modeling Group and collaborators

## Main principles:
*	All development should be shared in the group.
*	All changes to the code should be traceable.
*	We maintain just **one** code.

In order to share our developments and avoid that the code-branches drifts apart, we need to accept the risk that we sometimes push changes to developer branch that have bugs or create conflicts: **so rather update too often than too seldom.**  

## Maintaining the code
On a long term-basis we should keep only two branches:
1.	The master-branch
1.	The develop-branch

The master-branch (which should be stable) should be updated in the following cases:
*	For larger releases – e.g. update of the code fore delivery to met.no 
*	Bugs on the master should be fixed (as soon as possible) and the changes pushed to both the master and the developer branch.

The developer’s branch should be continuously updated:
*	All developments or bug fixes should be done on issue-branches, and those changes are pushed back to the develop-branch as soon as the issue is solved. After the issue is solved, the issue-branch should be discontinued.
*	Bug fixes on the master branch must also be pushed to the developer branch.

**If doubt when merging confer with Mostafa or Annette.**

## Quality control
### Current
When you have updated the code ask at least one person in the group to test (compile and do a small test-run) on a different machine.
### Planned for the future
We will have a standard test directory containing everything needed to perform a test test-run.
This directory will also contain scripts for standard validation/testing of the results.

## Version numbers
We will use version numbers, we have start with the a number 0.1.0 and go to 1.0.0 for the operational V4 upgrade. 

## Open issues on the GitHub
Each issue will be assigned to one person which is responible for solving and closing the issue.

## -----------------------------------------------------------------------------
## Usefull git commands for keeping your code up to date with the develop branch
## ----------------------------------------------------------------------


###  Download NERSC-HYCOM with git from Github:

> `$git clone https://github.com/nansencenter/NERSC-HYCOM-CICE.git`

****## How to track and update the 'develop' in GitHub from your local dir:****
* To show all local and remote branches:

` $ git branch -a  `

* To show only all remote branches:

`$ git branch -r       `

*  To track changes 'develop' branch in your local dir:

`$ git checkout --track origin/develop `
*  To switch to the 'develp' branch in your local dir: 

`$ git checkout develop   `      
        
* To update your local develop branch from the Github develop branch

`$ git pull origin develop              `


**>> Merging for example 'nesting' or 'Fram' in 'develop' branch and push it to the Github:**

1.  first switch to 'develop' 

`$ git checkout develop `

You can verify this by running '$git branch -a' and you will see a * pointing to 'develop' on the list 
2.  Now run merge command 

`$ git merge --no-ff nesting`


This will update 'develop' with all new changes in 'nesting' in you local dir:
If you there are some conflicts in some files, then you have 
to resolve them by opening these files one by one.

3. Now add and commit

`$ git add .`

`$ git commit -m ‘nesting is now merged to develop’`		# record changes in the new branch

To pull updates the corresponding 'develop' branch in Github:

` $ git push origin develop`



## # * General command Creating a new branch/edit and commit in your local dir and merge with develop branch:


* create a new branch


`$ git branch  featue1	`  
                                
*  switch to the new branch

`$ git checkout feature1	     `                            

**>> Now you edit files and make change, then add and commit:**

`$ git add .`

* record changes in the new branch

`$ git commit -m ‘add test’  `
                                
* compare with the develop branch

`$ git diff develop   `   
                           
*  switch back to the develop branch to be your active branch

`$git checkout develop `
                           
* merge  feature1

`$git merge --no-ff feature1     `                        
