# dse14-finding-venusian-volcanoes

This is the repository for all code related to the DSE of the Aerospace
Engineering education of the TU Delft in the academic year 2015-2016.

# How to work with git?

Mainly: Google it. But in short, there's two ways of doing this:
1. Download the software from git to work within a graphical environment
2. Use a terminal

I would advise you to use the terminal, because once you get used to it (which
is relatively quickly) you will work much faster in it. A really basic
introduction to the git workflow will be given below. If you want to see
visualisations or a better description, then see:
https://www.atlassian.com/git/tutorials/comparing-workflows/centralized-workflow

In git there are branches. To put it simply: a branch is a complete copy of all
code at a given point in time. All branches are given names. The main branch
(or, to use the existing analogy with trees) is sometimes called a trunk and
bears the name 'master' (or: development, but as we're a small team I'll just be
using the 'master' branch).

Lets say you want to implement a specific function in the project. You will then
consecutively take the following steps (the first term is the console command
that will achieve whatever you want to do):

1. `git checkout master` this will select the 'master' branch as the current
branch you're working on. This is necessary for the next command.
2. `git checkout -b my-new-branch-name` or two equivalent commands:
`git branch my-new-branch-name` followed by `git checkout my-new-branch-name`.
Both commands do the same thing: creating a copy of the branch you're currently
working on (which is the 'master' branch), naming it 'my-new-branch-name', and
then selecting it as your active branch.
3. Start coding for a bit, adding new stuff or editing old stuff.
4. `git add .` followed by `git commit -m 'this is a comment for my commit'`.
The first command (`git add`) will add files that you're going to commit (this
process is called 'staging'). The dot in the command `git add .` will indicate
that you want to add all files which were changed to the staging area. In case
you only want to add certain files then specify them by filename (e.g.
`git add ./folder/filename.py`). The second commit `git commit` will add them to
the local branch officially. Up until now you've been editing files, but they
are not yet associated with the branch you're working on. With `commit` you'll
do this. The added option (`-m 'comments'`) will associate a comment with you're
commit so its easy for yourself to find back earlier comments or for somebody
else to figure out in which steps you've been editing code.
5. Go to point (3) a couple of times. Do not be afraid to call `git add .` and
`git commit` often. If something goes wrong you can revert to a previous commit.
6. Once you're done perform `git push`. Usually git will give you an error
message containing the command you must be typing (with an `upstream`), retype
the given command to upload your changes to the repository (which is GitHub at
the moment).

You've now successfully created a branch, added code to it and uploaded all
changes to the repository. What remains to do is to make a `pull request`. This
can be done online on github. This will request the manager of the repository to
merge your changes with the main branch 'master'.

Several other important things you can do with git:

- `git merge another-branch` is a command you can use to merge the other branch
with the branch you're working on. For example. Lets say somebody made a pull
request and the manager merged the changes in. You still have an old version of
the master branch on your laptop and wish to retrieve the changes and merge them
into your branch, as you're working on something requiring the newly made
changes:
	1. `git checkout master`: selecting the master branch as the active branch
	2. `git pull`: retrieve/download the most recent version of the master branch
	3. `git checkout my-branch`: return to the branch you were working on that
	required the newly made changes to the master branch
	4. `git merge master`: merge all changes in the master branch into your local
	branch. Note: Sometime git will complain about 'merge conflicts'. This will mean
	that the automated system could not figure out which code to use from which
	branch. The program will automatically edit the files to include code from both
	branches. What you're supposed to do is to search through the entire project for
	the text `<<<<<` (indicating the start of a merge conflict), look at the code,
	decide which code needs to used and remove the other code. After that you can
	add, commit and push the resolved merge conflict.
- `git status` will give you a quick overview of all files that have been
changed and added.
- `git stash` will take all your current changes and put them onto a stack of
changes. All your files will revert to the last commit you made. Afterwards you
can do one of two things: `git stash pop`, which will put the changes which are
residing in the stack back into the files you're working on. This is useful when
you were accidentally working on the wrong branch (`git stash`,
`git checkout other-branch`, `git stash pop`, `git add .`,
`git commit -m 'now the code is on the correct branch'`
`git checkout my-branch`). Another possibility is `git stash drop`, which will
remove the changes which were put onto the stack. This is useful if you made
some changes in files and saved them, but wish to revert to a previous commit.
In that case you type `git stash` followed by `git stash drop` and all your
changes are gone. Not that this process cannot be reverted.

# SSH Keys, what are they?

Encryption keys ensure that the changes uploaded to GitHub are actually
coming from you and not from somebody else. In order to figure out how to
generate them:

https://confluence.atlassian.com/bitbucketserver/creating-ssh-keys-776639788.html

Once you've created an SSH key you can add it to your github account. Sometimes
you'll have to generate a new one when you're on a different network.

# IMPORTANT: In case of git mistakes

One can make mistakes in git. For example committing to the wrong branch,
pushing changes which shouldn't have been pushed to the remote branch (a remote
branch is the one residing on GitHub, a local branch is residing on your
computer, they're synchronized by calling `git push` or `git pull`, depending on
what you want to do). Usually they can be fixed by googling the problem, ending
up on stackoverflow and doing what people tell you to do. But if you're making
possibly disastrous changes, please ask somebody with knowledge of git to fix it
before trying anything!
On the other hand, if you're trying to fix a comment on a commit because your
grammar was incorrect, please don't bug somebody who has knowledge of git but
fix it yourself.
