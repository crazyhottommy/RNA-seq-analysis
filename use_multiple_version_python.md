
### Snakemake on my Mac

Recently, I want to test [`snakemake`](https://bitbucket.org/snakemake/snakemake/wiki/Documentation) for building up a pipeline. It requires python3, but I have python2.x installed.  
It seems [`virtualenv`](http://docs.python-guide.org/en/latest/dev/virtualenvs/) and [`pyenv`](https://github.com/yyuu/pyenv) are
the tools to manage multiple versions of python on the same machine. I tried `pyenv` on my mac OS.

install `pyenv`

```bash
brew update
brew install pyenv
```
add `if which pyenv > /dev/null; then eval "$(pyenv init -)"; fi` to my `.zshrc` file.

install python 

```bash
## what version of python are avaiable
pyenv install --list

## install python 3.5
pyenv install 3.5.1

## install python 2.7
pyenv install 2.7.11 

## set global 
pyenv global 3.5.1
pyenv versions
system
  2.7.11
* 3.5.1 (set by /Users/mtang1/.pyenv/version)

## if you want to use python 3.5.1 in the current folder
pyenv local 3.5.1
python --version
3.5.1

## reset to 2.7.11
pyenv local 2.7.11
python --version
2.7.11

```

install `snakemake`:

```bash
pyenv install miniconda3-3.19.0
pyenv global miniconda3-3.19.0
conda install -c bioconda snakemake
```

I got an error when evoke `snakemake` in the terminal:

>Traceback (most recent call last): File "/Users/mtang1/.pyenv/versions/miniconda3-3.19.0/bin/snakemake", line 6, in <module>   sys.exit(snakemake.main()) File "/Users/mtang1/.pyenv/versions/miniconda3-3.19.0/lib/python3.5/site-packages/snakemake/init.py", line   885, in main parser = get_argument_parser() File   "/Users/mtang1/.pyenv/versions/miniconda3-3.19.0/lib/python3.5/site-packages/snakemake/init.py", line 494, in get_argument_parser   const=available_cpu_count(), File "/Users/mtang1/.pyenv/versions/miniconda3-3.19.0/lib/python3.5/site-packages/snakemake/utils.py",   

I have to add `import multiprocessing` in the `utils.py` to avoid that.
[The issue](https://bitbucket.org/snakemake/snakemake/issues/324/nameerror-name-multiprocessing-is-not) was fixed in the new release.

### Snakemake on HPC

[Samir](https://github.com/dyndna) has installed `snakemake` on nautilus using `anaconda` and set up the virtual environment by `virtualenv`.
To activate it:

`source activate Rconda`

To inactivate it

`source deactivate`

### visualize `snakemake` workflow

```bash
## this will install dot 
brew install graphviz

snakemake --dag | dot | display
```
