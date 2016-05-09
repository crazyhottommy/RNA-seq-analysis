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


```
