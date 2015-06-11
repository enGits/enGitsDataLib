#!/usr/bin/sh
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# +                                                                    +
# + This file is part of enGitsDataLib.                                +
# +                                                                    +
# + Copyright 2015 enGits GmbH                                         +
# +                                                                    +
# + enGitsDataLib is free software: you can redistribute it and/or     +
# + modify it under the terms of the GNU Lesser General Public License +
# + as published by the Free Software Foundation, either version 3 of  +
# + the License, or (at your option) any later version.                +
# +                                                                    +
# + enGitsDataLib is distributed in the hope that it will be useful,   +
# + but WITHOUT ANY WARRANTY; without even the implied warranty of     +
# + MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the      +
# + GNU Lesser General Public License for more details.                +
# +                                                                    +
# + You should have received a copy of the GNU Lesser General Public   +
# + License along with enGitsDataLib.                                  +
# + If not, see <http://www.gnu.org/licenses/>.                        +
# +                                                                    +
# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
cd /tmp/edl/
doxygen
cp ./logo.png /srv/www/htdocs/edl
cd /tmp
rm -rf edl
