# Contributing to Bioconda

[Bioconda](https://bioconda.github.io/index.html) is a distribution of bioinformatics software for the Conda package manager. It is a community-driven project that provides a large collection of bioinformatics packages, making it easier for researchers to install and manage software dependencies.

[Instructions for contributing to Bioconda](https://bioconda.github.io/contributor/index.html) are available, but provide much more information than needed for researchers who intend to simply provide their package in Bioconda with no more care than to keep it updated. This tutorial will provide a simplified version of the instructions, focusing on the most important steps to create a package for a bioinformatics tool.

Before starting this tutorial, you should have a working version of your package, i.e., a [GitHub release](http://docs.github.com/en/repositories/releasing-projects-on-github/about-releases). If you don't, go to `https://github.com/<your_usename>/<your_package_name>/releases` and create a new release. This will create a tarball with the source code of your package, which is what we will use to build the bioconda package.

## Index

1. [Install conda](#1-install-conda)
2. [Setup your local project](#2-setup-your-local-project)
   * [Create a branch](#21-create-a-branch)
3. [Create a recipe](#3-create-a-recipe)
    * [meta.yaml](#31-metayaml)
        * [The build scripts](#the-build-scripts)
        * [About tests](#about-tests)
        * [About licensing](#about-licensing)
    * [build.sh](#32-buildsh)
4. [Add the files to your local copy of bioconda-recipes](#4-add-the-files-to-your-local-copy-of-bioconda-recipes)
5. [Test your recipe locally](#5-test-your-recipe-locally)
6. [Commit and push your changes](#6-commit-and-push-your-changes)
7. [Create a pull request and delete the branch](#7-create-a-pull-request-and-delete-the-branch)
8. [Keeping the package updated](#8-keeping-the-package-updated)

## 1. [Install conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

I suggest using [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/install), since Anaconda brings a lot to the table which is not needed for this tutorial (nor for most users).

Also, set the priority of the channels as follows:
```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```
This will make sure that the packages from Bioconda are preferred over the ones from conda-forge, which is important for the bioconda channel to work properly.

## 2. [Setup your local project](https://bioconda.github.io/contributor/setup.html)

First [create your own fork](https://github.com/bioconda/bioconda-recipes/fork) of the `bioconda-recipes` repository.

If you intend for the Scientific Computing ACF to manage your recipe, you can skip this step and use the [HCEMM fork](https://github.com/HCEMM/bioconda-recipes) instead.

Second, clone your fork to your local machine:

```bash
USERNAME=<your-username>
git clone https://github.com/${USERNAME}/bioconda-recipes
git clone https://github.com/HCEMM/bioconda-recipes         # If you want to use the HCEMM fork
cd bioconda-recipes
git remote add upstream https://github.com/bioconda/bioconda-recipes        # add the main bioconda-recipes repo as an upstream remote - this will simply inform git that this is the original repo, and make it easier to keep your fork up to date
```

### 2.1. Create a branch

Create a new branch for your changes. This is a good practice, as it allows you to work on your changes without affecting the main branch.

```bash
# Make sure our master is up to date with Bioconda
git checkout master
git pull upstream master
git push origin master

# Create and checkout a new branch for our work
git checkout -b <add-my-recipe>
```

## 3. [Create a recipe](https://bioconda.github.io/contributor/workflow.html#id3)

To keep it simple, we provide here an example that will be appropriate for the most common cases.

### 3.1. meta.yaml

The [`meta.yaml`](https://docs.conda.io/projects/conda-build/en/stable/resources/define-metadata.html) file defines the metadata for the package, and is the essential piece of the recipe. The metadata includes the package name, version, dependencies, and other information.

Here follows an example `meta.yaml` file:

```yaml
{% set username = "<username>" }
{% set name = "<package_name>" %}
{% set version = "<0.1.0>" %}
{% set sha256 = "<SHA256>" %}

package:
  name: {{ name|lower }}    # set the name of the package (to lowercase)
  version: {{ version }}    # set the version of the package

source:
  url: https://github.com/{{ username }}/{{ name }}/archive/{{ version }}.tar.gz    # URL for the source code
  sha256: {{ sha256 }}      # SHA256 hash of the source code tarball - explained below

build:
  noarch: generic     # set to generic if the package is architecture-independent, otherwise set to the architecture (e.g. linux-64, osx-64, win-64)
  number: 0           # set the build number - if the build fails publicly, increment this number in shame and try again
  run_exports:        # set the run exports - no idea what this does, but is necessary ¯\_(ツ)_/¯
    - { { pin_subpackage(name, max_pin="x.x") } }
  script: >           # the build steps, explained below
    mkdir -p $PREFIX/bin && 
    cp <myscript> $PREFIX/bin &&
    chmod +x $PREFIX/bin/<myscript> &&
    ln -s $PREFIX/bin/<myscript> $PREFIX/bin/{{ name }}

requirements:
  run:                # required dependencies for the package, see: https://docs.conda.io/projects/conda-build/en/stable/resources/define-metadata.html#run
    - <dependency1>
    - <dependency2>

test:                 # tests to run after the package is built. See below.
  commands:
    - {{ name }} -v

about:
  home: https://github.com/{{ username }}/{{ name }}      # the homepage of the package. If you don't have one, set it to the GitHub page, with a cool README
  license: BSD-3-Clause       # the license of the package - see below
  license_family: BSD
  license_file: LICENSE       # name of the License file - see below
  summary: '<This tool does something useful>'      # a short summary of the package
  description: |          # a less short summary of the package          
    <This tool does something useful for bioinformatics.>
    <And I am proud of it.>
  doc_url: https://github.com/{{ username }}/{{ name }}/wiki      # link to the documentation of the package. Link to the README.md if you feeling lazy
  dev_url: https://github.com/{{ username }}/{{ name }}           # link to the development page of the package. Usually the GitHub page

extra:
  recipe-maintainers:
    - {{ username }}      # the maintainer of the recipe. You!, or if you are using the HCEMM fork, set to "iquasere"
```

Values in `<>` should be replaced with the actual values for your package:
- `username`: Your GitHub username. This is used to create the URL for the source code.
- `package_name`: The name of your package. This should be the same as the name of the repository on GitHub.
- `version`: The version of your package. This should be the same as the version tag on GitHub.
- `sha256`: The SHA256 hash of the tarball. This is used to verify the integrity of the package. You can obtain this value by running 
```bash
PACKAGE=<your-package-name>
VERSION=<your-package-version>
wget -O- https://github.com/${USERNAME}/${PACKAGE}/archive/refs/tags/${VERSION}.tar.gz | shasum -a 256
```

#### The build scripts

The `build.sh` and `bld.bat` files are used to build the package for Linux and Windows, respectively. These files contain the commands that will be run to build the package, and are better explained [here](https://docs.conda.io/projects/conda-build/en/stable/resources/build-scripts.html). These build commands can also be included in the `meta.yaml` file if simple enough, but otherwise keep them to the build scripts.

These steps follow the reallocation of the files from inside the GitHub repo to the right places in the conda environment. There are many implicit things here, e.g., your tarball has been untared, and you are at the root of your GitHub repo. The `$PREFIX` variable here is essential, and denotes the root of the conda environment.

For your tool to be available from anywhere, its executable must be on the `bin` folder. Here, `mkdir -p $PREFIX/bin` assures that folder exists, and `cp <myscript> $PREFIX/bin` copies the main script to that folder. `chmod +x $PREFIX/bin/<myscript>` makes the script executable, and `ln -s $PREFIX/bin/<myscript> $PREFIX/bin/{{ name }}` creates a symbolic link to the script in the `bin` folder so, for example, `myamazingtool.py` can be called as `myamazingtool`. It will also be able to be called as `myamazingtool.py`, because that script is in `$PREFIX/bin`.
```yaml
script: >
  mkdir -p $PREFIX/bin && 
  cp <myscript> $PREFIX/bin &&
  chmod +x $PREFIX/bin/<myscript> &&
  ln -s $PREFIX/bin/<myscript> $PREFIX/bin/{{ name }}
```
If you have a lot of files to copy, including support tables or other scripts, put them on `$PREFIX/share`. `bin` should be kept for the main files only.

#### About tests

Here I put a simple test that checks if the main script is callable, and requests the version from it (with the `-v` flag). This is a good test to check if the package is accessible, and it can import its dependencies, i.e., whatever dependencies are imported in the main script.

```yaml
test:
  commands:
    - {{ name }} -v
```
If you want to add additional tests (you should), don't put them here, as this will strain Bioconda's CI. Instead, [use GitHub Actions](https://github.com/HCEMM/sc_tutorials/blob/main/github_tutorial/README.md#3-github-actions-cicd), or some other CI service.

#### About licensing

To choose the license that best suits your needs, see [GitHub's guidelines](https://docs.github.com/en/repositories/managing-your-repositorys-settings-and-features/customizing-your-repository/licensing-a-repository
) and [ours](https://github.com/HCEMM/sc_tutorials/blob/main/github_tutorial/README.md#5-licensing-and-citation).

```yaml
license_file: LICENSE
```
The parameter for `license_file` is the exact name of the license file in the repository. If the file is at the root and is called `LICENSE`, this is the right value. Adjust otherwise.

### 3.2. build.sh

If your package requires more than a few commands to set it up - you have many helper files, or many scripts changing places -, the recipe is better implemented using a `build.sh` file. Here follows the same example as above, but using a `build.sh` file instead of the `script` parameter in the `meta.yaml` file.

In the `meta.yaml`, remove the `script` parameter from the `build` section. The `build.sh` file will be used instead, Bioconda will recognize it automatically. The `build.sh` file then becomes:

```bash
PACKAGE_NAME = <package_name>
MYSCRIPT = <myscript>
mkdir -p $PREFIX/bin && 
cp $MYSCRIPT $PREFIX/bin &&
chmod +x $PREFIX/bin/<myscript> &&
ln -s $PREFIX/bin/$MYSCRIPT $PREFIX/bin/$NAME
```
Again, consider you are at the root of your repo, and `$PREFIX` is the root of the conda environment.

## 4. Add the files to your local copy of bioconda-recipes

Go to your local copy of `bioconda-recipes` repository, create a folder for your package, and put the files created previously inside there. The folder name should be the same as the package name, in lowercase. This is important, as Bioconda uses this folder to find the recipe.
```bash
cd bioconda-recipes
mkdir recipes/$PACKAGE_NAME
cp meta.yaml build.sh recipes/$PACKAGE_NAME
```

## 5. [Test your recipe locally](https://bioconda.github.io/contributor/building-locally.html)

Bioconda runs a complex workflow to avoid including disfunctional recipes. This workflow looks for your tarball (specified with the `url` parameter in the `meta.yaml` file) and builds the package from it. It also runs the tests specified, but again, this is better relegated to GitHub Actions or some other CI service, to not strain Bioconda's CI.

In this section, we will test the recipe locally, using the `conda-build` command. This command will build the package from the recipe, and run the tests specified in the `meta.yaml` file. This is a good way to check if the recipe is working properly before submitting it to Bioconda, saving time for you and resources from Bioconda.

```bash
mamba create -n bioconda -c conda-forge -c bioconda bioconda-utils

conda activate bioconda

# linting assures the sections of the meta.yaml file are correct
bioconda-utils lint --git-range master

# build a testing environment, similar to Bioconda's CI
bioconda-utils build --docker --mulled-test --git-range master
```

## 6. [Commit and push your changes](https://bioconda.github.io/contributor/workflow.html#id4)

This will put your changes on your fork of the `bioconda-recipes` repository, remotely (at GitHub).

```bash
git add *
git commit -m "Add <package name> recipe"   # remember to change to your package name!
git push
```

## 7. Create a pull request and delete the branch

[Creating the pull request](https://bioconda.github.io/contributor/workflow.html#id6) will trigger the build process of Bioconda. If everything goes alright, your tool will be available for install with the `conda install` command in a few hours.

[Delete the branch](https://bioconda.github.io/contributor/workflow.html#id7) to keep your fork clean. You wont need this branch anymore.

## 8. Keeping the package updated

This is the fun part. When you launch new versions on Github, a pull request will be automatically created for `bioconda-recipes`! After being accepted, your tool will be available through `conda install` with the new code and version.

This is if the build script remains the same. If, by virtue of your tool evolution, it requires extra steps on the build script, or extra dependencies. these can be added manually in the pull request.
