FROM gitpod/workspace-full

USER gitpod

RUN pyenv install 2.7.15 && pyenv global 2.7.15 3.6.6 && pip install ipython jupyter

RUN pip install sympy