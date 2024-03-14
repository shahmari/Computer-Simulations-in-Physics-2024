# Getting Started with Git

## Introduction

Git is a powerful version control system widely used by developers to manage source code and collaborate on projects. It allows you to track changes, revert to previous versions, and work on different features simultaneously. This tutorial will guide you through the basics of using Git.

## Prerequisites

Before getting started, ensure that Git is installed on your system. You can download and install Git from the official website: [Git - Downloads](https://git-scm.com/downloads)

## Setting Up Git

Once Git is installed, open a terminal or command prompt and configure Git with your name and email address. This information will be attached to your commits.

```bash
git config --global user.name "Your Name"
git config --global user.email "your.email@example.com"
```

## Creating a New Repository

To start using Git, you need to create a new repository. Navigate to the directory where you want to store your project and initialize a new Git repository.

```bash
cd path/to/your/project
git init
```

This command creates a new Git repository in the current directory.

## Adding Files

After initializing the repository, you can start adding files to it. Create a new file or copy existing files into the project directory. To add files to the staging area (preparing them for the next commit), use the `git add` command.

```bash
git add filename
```

Replace `filename` with the name of the file you want to add. You can also use `git add .` to add all files in the current directory.

## Committing Changes

Once you have added the necessary files to the staging area, you can commit them to the repository. Commits are snapshots of your project at a specific point in time.

```bash
git commit -m "Commit message"
```

Replace `"Commit message"` with a brief description of the changes you made in this commit.

## Pushing Changes to a Remote Repository

After committing your changes to your local repository, you may want to share them with others or back them up on a remote server (like GitHub, GitLab, or Bitbucket). The `git push` command is used to upload local repository content to a remote repository.

```bash
git push origin branchname
```

Replace `branchname` with the name of the branch you want to push to the remote repository. `origin` refers to the default name of the remote repository, but you can replace it with any custom remote name if needed.

## Pulling Changes from a Remote Repository

When working in a collaborative environment, it's essential to keep your local repository up-to-date with changes made by others. The `git pull` command is used to fetch and integrate changes from a remote repository into your local repository.

```bash
git pull origin branchname
```

This command will fetch changes from the `origin` remote repository and merge them into the current branch in your local repository.

## Cloning a Repository

If you want to start working on an existing project hosted on a remote repository, you can clone it to create a local copy on your machine. Use the `git clone` command followed by the URL of the remote repository.

## Viewing Commit History

To view the commit history of a repository, you can use the `git log` command. This command displays a list of commits in reverse chronological order.

```bash
git log
```

You can customize the output of `git log` using various options like `--oneline`, `--graph`, and more. Refer to the Git documentation for more details.

```bash
git clone repository_url
```

This command will create a new directory containing a copy of the remote repository.

## Viewing Repository Status

At any time, you can check the status of your repository to see which files have been modified, staged, or committed.

```bash
git status
```

This command displays the current status of your repository.

## Branching and Merging

Git allows you to work on different features or versions of your project simultaneously using branches. To create a new branch, use the `git branch` command.

```bash
git branch branchname
```

Replace `branchname` with the name of your new branch. You can switch between branches using `git checkout`.

```bash
git checkout branchname
```

Once you have made changes in a branch and want to merge those changes back into the main branch (usually `master`), use the `git merge` command.

```bash
git checkout master
git merge branchname
```

This will merge the changes from `branchname` into the `master` branch.

## Ignoring Files

Sometimes, there are files or directories in your project that you don't want Git to track, such as temporary files, build artifacts, or sensitive information. You can specify patterns for Git to ignore using a special file called `.gitignore`.

Create a file named `.gitignore` in the root directory of your repository and list the file patterns you want Git to ignore.

Example `.gitignore` file:
```plaintext
*.log
build/
credentials.txt
```

This file tells Git to ignore all `.log` files, the `build` directory, and `credentials.txt`.

## Conclusion

This tutorial covers the basic commands and concepts of Git to help you get started with version control. As you become more familiar with Git, you can explore more advanced features and workflows to streamline your development process.