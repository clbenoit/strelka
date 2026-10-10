#
# Strelka - Small Variant Caller
# Copyright (c) 2009-2018 Illumina, Inc.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#

# gcc 14 -O3 false positive in boost random headers (-Wmaybe-uninitialized on the
# boost::binomial_distribution ctor): disable that single warning for this
# simulation-only translation unit (mechanism provided by cxxLibrary.cmake)
set(starling_sim_test_COMPILE_FLAGS "-Wno-maybe-uninitialized")

include(${THIS_CXX_LIBRARY_CMAKE})
