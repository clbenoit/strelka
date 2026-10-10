#!/usr/bin/env bash

#
# This documents the conversion of full boost to the subset contained in this project
#
# note this is less precise than what you would expect from bcp... the goal is to leave most scientifically practical
# header-only libraries in place, and remove large or infrequently used cases
#

set -o nounset
set -o xtrace

rel2abs() {
  pwd -P $1
}

thisDir=$(rel2abs $(dirname $0))

mkdir -p output
cd output

# boost version (argument 1, default 1.92.0): the full tarball must be in this script's directory as $boost_name.tar.bz2
boost_name=${1:-boost_1_92_0}
output_name=${boost_name}_subset

boost_tarball=$thisDir/$boost_name.tar.bz2

if ! [ -f $boost_tarball ]; then
    echo "Can't find input boost tarball '$boost_tarball'"
    exit 1
fi

tar -xjf $boost_tarball
mv $boost_name $output_name

for ddir in doc more status; do
    rm -rf $output_name/$ddir
done

# remove unused headers for larger sub-libraries:
for f in asio geometry gil graph polygon python signals signals2; do
    rm -rf $output_name/boost/$f*
done

# remove unused libs:
# compiled libs kept in full (built by bjam per src/cmake/boost.cmake), plus wave kept for parity with the 1.58 subset:
libs_full="config date_time detail serialization timer chrono filesystem program_options system test wave"
# boost 1.92's b2 resolves project references declared in each lib's build.jam
# (e.g. date_time depends on /boost/algorithm) and fails if libs/<dep> is missing:
# keep the transitive header-only dependency closure, trimmed of test/example/doc
# (the headers themselves ship in the top-level boost/ tree):
libs_headers="algorithm align any array assert atomic bind compat concept_check container container_hash conversion core describe endian exception filesystem function function_types fusion headers integer intrusive iterator lexical_cast move mp11 mpl optional pool predef preprocessor proto range ratio regex scope serialization smart_ptr spirit system test thread throw_exception timer tokenizer tuple type_index type_traits typeof unordered utility variant variant2 winapi io phoenix"
# io and phoenix have no build.jam in 1.92 (plain header-only dirs) but the dirs must
# exist for b2 project resolution
# the numeric_conversion project is nested: Jamroot maps it to libs/numeric/conversion;
# keep that subtree, trim the other nested libs (ublas, interval, odeint, ...):
(
cd  $output_name/libs/numeric && ls | grep -v -E "^(conversion)$" | xargs rm -rf
)
rm -rf $output_name/libs/numeric/conversion/test $output_name/libs/numeric/conversion/example
(
cd  $output_name/libs
ls | grep -v -E "^($(echo $libs_full | tr ' ' '|')|numeric)$"  | grep -v -E "^($(echo $libs_headers | tr ' ' '|'))$"  | xargs rm -rf
)

# trim test/example of the header-only dependency libs (not of the compiled seeds):
trim_list=""
for f in $libs_headers; do
    case " $libs_full " in *" $f "*) ;; *) trim_list="$trim_list $f" ;; esac
done
for f in $trim_list; do
    rm -rf $output_name/libs/$f/test $output_name/libs/$f/example
done

# remove unused tools:
# boost_install is required since boost 1.77+: the b2 `install` rules live in tools/boost_install/boost-install.jam
(
cd  $output_name/tools
ls | grep -v -e "^\(build\|boost_install\|inspect\)$"  | xargs rm -rf
)

# remove docs:
(
cd $output_name
find . -name doc -type d -print | xargs rm -rf
)

# tarball up:

tar -c $output_name -f - | bzip2 -c -9 >| $output_name.tar.bz2
