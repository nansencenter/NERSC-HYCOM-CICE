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
