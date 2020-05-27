"""
Workaround for https://github.com/bitprophet/releases/issues/65

This has to be in its own file, as for some reason classes in conf.py cannot
be pickled.
"""

class release_uri:
    def __init__(self, releases_github_path):
        self._path = releases_github_path

    def __mod__(self, release):
        if release[0].isdigit():
            release = "v" + release
        return 'https://github.com/%s/tree/%s' % (self._path, release)
