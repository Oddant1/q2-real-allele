# flake8: noqa
# ----------------------------------------------------------------------------
# Copyright (c) 2026, Oddant1.
#
# No license has been chosen for this software yet.
#
# See the LICENSE file distributed with this software.
# ----------------------------------------------------------------------------

try:
    from ._version import __version__
except ModuleNotFoundError:
    __version__ = '0.0.0+notfound'
