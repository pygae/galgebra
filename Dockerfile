FROM gitpod/workspace-full

USER gitpod

RUN pyenv install 2.7.15 && pyenv global 3.6.6 2.7.15

RUN pip2 install ipython jupyter && pip3 install ipython jupyter

RUN pip2 install sympy && pip3 install sympy

RUN pip2 install pytest nbval pytest-cov && pip3 install pytest nbval pytest-cov