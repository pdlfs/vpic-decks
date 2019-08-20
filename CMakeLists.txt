#
# Copyright (c) 2019 Carnegie Mellon University,
# Copyright (c) 2019 Triad National Security, LLC, as operator of
#     Los Alamos National Laboratory.
#
# All rights reserved.
#
# Use of this source code is governed by a BSD-style license that can be
# found in the LICENSE file. See the AUTHORS file for names of contributors.
#

#
# CMakeLists.txt  deltafs-umbrella scripts
# 20-Aug-2019  chuck@ece.cmu.edu
#

# command line config:
#
# -DCMAKE_INSTALL_PREFIX=/tmp/delta       # where to install
#                                         # defaults to /usr/local
#

cmake_minimum_required (VERSION 3.0)
project (vpic-decks C CXX)

# decks to install
set (install-decks forcefree-tracing)

#
# "make install" rules
#
install (DIRECTORY ${install-decks} DESTINATION decks)
