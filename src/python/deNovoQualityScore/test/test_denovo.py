#!/usr/bin/env python3
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

import gzip
import itertools
import os
import re
import subprocess as sp
import sys
import warnings

import cyvcf2
import numpy as np
import pytest

from numpy.testing import assert_allclose

import denovo


python_exec = sys.executable
denovo_script = './denovo.py'

dev_null = open(os.devnull, 'wb')

strelka_truth_set = [
    ('chr20', 54192843, 'G', ['A']),
    ('chr20', 54455224, 'G', ['C']),
    ('chr20', 55297452, 'A', ['T']),
    ('chr20', 60000000, 'A', ['AACACACACAC']),  # simulated
    ('chr20', 60000001, 'ATC', ['A', 'ATCTC'])  # simulated
]

# manually inspected
manta_truth_set = [
    ('chr22', 10686252, 'BND'),
    ('chr22', 12188286, 'BND'),
    ('chr22', 12397725, 'DEL'),
    ('chr22', 22266782, 'BND'),
    ('chr22', 23510376, 'BND'),
    ('chr22', 30597450, 'BND'),
    ('chr22', 31355861, 'DEL'),
    ('chr22', 36395692, 'INV'),
    ('chr22', 39107707, 'DEL')
]


def get_test_data_path(file_name, dir_name='data', space_name=''):
    """Compute full path to a test file"""

    script_dir = os.path.abspath(os.path.dirname(__file__))
    rel_path = os.path.join(space_name, dir_name, file_name)
    path = os.path.join(script_dir, rel_path)

    if not os.path.exists(path):
        msg = "File '{0}' does not exist.".format(path)
        raise IOError(msg)

    return path


@pytest.fixture(scope='session')
def datasets():

    files = {
        'vcf_in': 'denovo-chr1-200-snv.vcf.gz',
        'vcf_in_phased': 'denovo-chr1-200-snv-phased.vcf.gz',
        'vcf_spw': 'variants-strelka-slice.vcf.gz',
        'vcf_spw_5samples': 'variants-strelka-five-fake-samples.vcf.gz',
        'vcf_spw_chr22': 'spw312-pedigree-chr22.vcf.gz',
        'vcf_spw_sv_chr22': 'pedigree-sv-chr22.vcf.gz',
        'vcf_spw_chrXY': 'spw-allosomes-male-proband.vcf.gz',
        'vcf_spw_pedphase': 'spw-791.vcf.gz',
        'vcf_spw542': 'spw542-cases.vcf.gz',
        'vcf_spw1005': 'spw1005-prior.vcf.gz',
        'invalid_vcf': 'denovo-chr1-200-snv.ped',
        'variant_examples': 'variant-examples.vcf.gz',
        'dng_ref': 'dng-auto-ref.txt.gz',
        'mutation_AT': 'mutation-AT.vcf',
        'mutation_AT_minimal': 'mutation-AT-minimal.vcf',
        'missing_PL': 'mutation-AT-missing-PL.vcf',
        'missing_PL': 'mutation-AT-missing-PL.vcf',
        'missing_DP': 'mutation-AT-missing-DP.vcf',
        'par_vcf': 'spw588-par.vcf.gz',
        'par_bed': 'PARv5.bed'
    }

    paths = {name: get_test_data_path(path) for name, path in files.items()}

    return paths


@pytest.fixture(scope='function')  # new tempdir for each test
def tempdir(tmpdir_factory):

    # set-up
    tmp_dir = tmpdir_factory.mktemp('tmp')

    # passing
    yield tmp_dir

    # tear-down
    tmp_dir.remove()


dev_null = open(os.devnull, 'wb')

strelka_truth_set = [
    ('chr20', 54192843, 'G', ['A']),
    ('chr20', 54455224, 'G', ['C']),
    ('chr20', 55297452, 'A', ['T']),
    ('chr20', 60000000, 'A', ['AACACACACAC']),  # simulated
    ('chr20', 60000001, 'ATC', ['A', 'ATCTC'])  # simulated
]

# manually inspected
manta_truth_set = [
    ('chr22', 10686252, 'BND'),
    ('chr22', 12188286, 'BND'),
    ('chr22', 12397725, 'DEL'),
    ('chr22', 22266782, 'BND'),
    ('chr22', 23510376, 'BND'),
    ('chr22', 30597450, 'BND'),
    ('chr22', 31355861, 'DEL'),
    ('chr22', 36395692, 'INV'),
    ('chr22', 39107707, 'DEL')
]


# commonly used functions

#
# VCF access helpers (cyvcf2):
# the test suite previously relied on PyVCF (python2-only). The helpers below
# provide the small subset of that API used here on top of cyvcf2, which is
# already a runtime dependency of denovo.py.
#
def get_denovo_program(vcf_reader):
    """value of the ##denovo_program header line (or None)"""
    for line in vcf_reader.raw_header.splitlines():
        if line.startswith('##denovo_program='):
            return line.split('=', 1)[1]
    return None


def get_formats(vcf_reader):
    """set of FORMAT field names declared in the VCF header"""
    fields = set()
    for line in vcf_reader.raw_header.splitlines():
        if line.startswith('##FORMAT='):
            for m in re.finditer(r'<ID=([A-Za-z0-9_]+)', line):
                fields.add(m.group(1))
    return fields


def is_snp(variant):
    """single nucleotide variant (PyVCF 'is_snp' semantics, simplified)"""
    return len(variant.REF) == 1 and all(len(alt) == 1 for alt in variant.ALT)


def is_indel(variant):
    """insertion/deletion (PyVCF 'is_indel' semantics, simplified)"""
    if 'SVTYPE' in variant.INFO:
        return False
    return len(variant.REF) > 1 or any(len(alt) != len(variant.REF) for alt in variant.ALT)


def is_sv(variant):
    """structural variant (PyVCF 'is_sv': SVTYPE INFO field present)"""
    return variant.INFO.get('SVTYPE') is not None


def var_type(variant):
    """'snp' / 'indel' / 'sv' / 'unknown' (PyVCF 'var_type' semantics)"""
    if is_snp(variant):
        return 'snp'
    elif is_indel(variant):
        return 'indel'
    elif is_sv(variant):
        return 'sv'
    else:
        return 'unknown'


def approx(a, b, rel_tol=1e-4, abs_tol=1e-6):
    """numeric equality tolerant to cyvcf2's float32 storage of FORMAT
    fields (PyVCF reported float64 values)"""
    fa, fb = float(a), float(b)
    return abs(fa - fb) <= max(rel_tol * max(abs(fa), abs(fb)), abs_tol)


def info_eq(info1, info2):
    """cyvcf2 INFO objects compare by identity: compare by content"""
    if info1 is None or info2 is None:
        return info1 is None and info2 is None
    return dict(info1) == dict(info2)


def variants_match(v1, v2):
    """field-by-field comparison replacing PyVCF's variant equality"""
    return all(getattr(v1, f) == getattr(v2, f)
               for f in ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'FILTER')) \
        and info_eq(v1.INFO, v2.INFO)


def has_field(variant, field):
    """True if the FORMAT field is present in the variant record"""
    return field in variant.FORMAT


def sample_field(variant, sample_idx, field):
    """single-valued FORMAT field for a sample; None if absent or missing.

    cyvcf2 reports missing values as NaN (numeric fields) or '.'
    (string fields), and raises KeyError for a field not in FORMAT;
    PyVCF reported them as None.
    """
    if not has_field(variant, field):
        return None
    value = variant.format(field)[sample_idx]
    if isinstance(value, str):
        # string field: value is a scalar (not a length-1 array)
        return value if value != '.' else None
    value = value[0]
    if value != value:  # NaN -> missing
        return None
    if isinstance(value, (int, float, np.integer, np.floating)):
        # cyvcf2 stores numeric FORMAT fields as float32: promote to
        # python float for py2-compatible comparison semantics
        return float(value)
    return value


def check_denovo_header(variant_reader):

    # test if output VCF contains the denovo program name and version
    reference_value = '{} {}'.format(denovo._name, denovo._version)
    assert get_denovo_program(variant_reader) == reference_value

    assert denovo._field_name in get_formats(variant_reader)

    assert denovo._field_name in get_formats(variant_reader)


# cmd interface

def test_cmd_missing_cmd():
    """Test cmd for non-existing denovo script"""

    cmd = ['tool_that_does_not_exist']

    with pytest.raises(OSError):
        sp.check_call(cmd, stdout=dev_null, stderr=dev_null)


def test_cmd_missing_vcf():
    """Test cmd for unspecified VCF file"""

    cmd = [python_exec, denovo_script]

    with pytest.raises(sp.CalledProcessError) as e:
        sp.check_call(cmd, stdout=dev_null, stderr=dev_null)
    assert e.value.returncode == 2


def test_premature_ended_pipe(datasets):
    """Test premature closure of output pipe"""

    n = 20
    cmd = [python_exec, denovo_script, '--proband', 'PROBAND', '--mother',
           'PARENT1', '--father', 'PARENT2', '--model', 'dummy', datasets['vcf_in']]
    cmd = ' '.join(cmd) + ' | head -n ' + str(n)

    p = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)

    # 'n' line of output
    assert len(p.stdout.readlines()) == n
    # no 'broken pipe' or no lines in stderr
    std_err = p.stderr.readlines()
    assert len(std_err) == 0 or not std_err.count('broken pipe')


def test_cmd_ids_full(datasets, tempdir):
    """Test cmd dummy model"""

    vcf_out_path = tempdir.join('out.vcf').strpath
    cmd = [python_exec, denovo_script, '--proband', 'PROBAND', '--mother',
           'PARENT1', '--father', 'PARENT2', '--model', 'dummy', datasets['vcf_in']]

    with open(vcf_out_path, 'w') as f:
        sp.check_call(cmd, stdout=f)

    vr = cyvcf2.VCF(vcf_out_path)

    check_denovo_header(vr)

    for variant in vr:
        assert 'DQ' in variant.FORMAT
        assert sample_field(variant, 0, 'DQ') == 0.0


def test_cmd_vcf_output_options(datasets, tempdir):
    """Test cmd interface with all VCF output options"""

    cmd = [python_exec, denovo_script,
           '--proband', 'Child', '--mother', 'Mother', '--father', 'Father',
           datasets['vcf_spw']]

    # header only
    cmd1 = cmd + ['--vcf-output', 'header']
    vcf_out_path1 = tempdir.join('out1.vcf').strpath
    with open(vcf_out_path1, 'w') as f:
        sp.check_call(cmd1, stdout=f)
    vr1 = cyvcf2.VCF(vcf_out_path1)

    check_denovo_header(vr1)
    assert len(list(vr1)) == 0

    # body only
    cmd2 = cmd + ['--vcf-output', 'body']
    vcf_out_path2 = tempdir.join('out2.vcf').strpath
    with open(vcf_out_path2, 'w') as f:
        sp.check_call(cmd2, stdout=f)
    # body-only output has no header, so it cannot be opened by cyvcf2:
    # check the header lines directly.
    with open(vcf_out_path2) as f:
        body_lines = [line for line in f if line.strip()]
        assert not any(line.startswith('##FORMAT=') for line in body_lines)
        assert not any(line.startswith('##INFO=') for line in body_lines)
        assert not any(line.startswith(('##contig=', '##sequence=')) for line in body_lines)
        assert len([line for line in body_lines if not line.startswith('#')]) > 0

    # full VCF
    cmd3 = cmd + ['--vcf-output', 'full']
    vcf_out_path3 = tempdir.join('out3.vcf').strpath
    with open(vcf_out_path3, 'w') as f:
        sp.check_call(cmd3, stdout=f)
    vr3 = cyvcf2.VCF(vcf_out_path3)

    check_denovo_header(vr3)
    assert len(list(vr3)) > 0

    # invalid option
    cmd4 = cmd + ['--vcf-output', 'definitely-not-an-option']
    with pytest.raises(sp.CalledProcessError) as e:
        sp.check_call(cmd4, stdout=dev_null, stderr=dev_null)
    assert e.value.returncode == 2


def test_cmd_ids_spw_ref(datasets, tempdir):
    """Test cmd interface against reference values"""

    cmd1 = [python_exec, denovo_script, '--proband', 'Child', '--mother', 'Mother',
            '--father', 'Father', datasets['vcf_spw']]
    # order of input arguments irrelevant
    cmd2 = [python_exec, denovo_script, '--proband', 'Child', '--mother', 'Mother',
            datasets['vcf_spw'], '--father', 'Father']

    vcf_out_path1 = tempdir.join('out1.vcf').strpath
    with open(vcf_out_path1, 'w') as f:
        sp.check_call(cmd1, stdout=f)

    vcf_out_path2 = tempdir.join('out1.vcf').strpath
    with open(vcf_out_path2, 'w') as f:
        sp.check_call(cmd2, stdout=f)

    vr1 = cyvcf2.VCF(vcf_out_path1)
    vr2 = cyvcf2.VCF(vcf_out_path2)

    check_denovo_header(vr1)
    check_denovo_header(vr2)

    n_hit = 0
    n_called = {'indel': 0, 'snp': 0}
    for variant, variant2 in zip(vr1, vr2):
        assert variants_match(variant, variant2)
        assert has_field(variant, 'DQ')
        v = (variant.CHROM, variant.POS, variant.REF, variant.ALT)
        dq = sample_field(variant, 2, 'DQ')
        if v in strelka_truth_set:
            assert dq >= 20.0
            n_hit += 1
        if dq is not None and dq >= 15.0:
            n_called[var_type(variant)] += 1
        # check number of decimals
        if dq is not None:
            assert approx(dq, round(dq, denovo._default_dq_decimals))
            assert dq <= denovo._max_score

    assert n_called['snp'] >= 4
    assert n_called['indel'] >= 2
    # check if all truth set variants were called: full recall
    assert n_hit == len(strelka_truth_set)


def test_cmd_parallel_spw_ref(datasets, tempdir):
    """Test parallel cli against reference values"""

    vcf_out_path1 = tempdir.join('out1.vcf').strpath
    cmd1 = [denovo_script, '--parallel',
            '--proband', 'Child',
            '--mother', 'Mother',
            '--father', 'Father',
            datasets['vcf_spw'], vcf_out_path1]

    sp.check_call(cmd1)

    # run the same analysis a second time on the previous output should not change the results
    vcf_out_path2 = tempdir.join('out2.vcf').strpath
    cmd2 = [denovo_script, '--parallel',
            '--proband', 'Child',
            '--mother', 'Mother',
            '--father', 'Father',
            vcf_out_path1, vcf_out_path2]

    #cmd2 = [denovo_script_parallel, 'Child', 'Mother', 'Father', vcf_out_path1, vcf_out_path2]

    sp.check_call(cmd2)

    vcf_out_path2 = vcf_out_path1

    vr1 = cyvcf2.VCF(vcf_out_path1)
    vr2 = cyvcf2.VCF(vcf_out_path2)

    check_denovo_header(vr1)
    check_denovo_header(vr2)

    n_hit = 0
    n_called = {'indel': 0, 'snp': 0}
    for variant, variant2 in zip(vr1, vr2):
        assert variants_match(variant, variant2)
        assert has_field(variant, 'DQ')
        v = (variant.CHROM, variant.POS, variant.REF, variant.ALT)
        dq = sample_field(variant, 2, 'DQ')
        if v in strelka_truth_set:
            assert dq >= 20.0
            n_hit += 1
        if dq is not None and dq >= 15.0:
            n_called[var_type(variant)] += 1
        # check number of decimals
        if dq is not None:
            assert approx(dq, round(dq, denovo._default_dq_decimals))
            assert dq <= denovo._max_score

    assert n_called['snp'] >= 4
    assert n_called['indel'] >= 2
    # check if all truth set variants were called: full recall
    assert n_hit == len(strelka_truth_set)


def test_cmd_single_parallel_spw_comparison(datasets, tempdir):
    """Test parallel cli against reference values"""

    # call 'denovo.py' directly, output to stdout
    vcf_out_path1 = tempdir.join('out1.vcf').strpath
    cmd1 = [denovo_script,
            '--proband', 'Child',
            '--mother', 'Mother',
            '--father', 'Father',
            datasets['vcf_spw']]

    with open(vcf_out_path1, 'w') as f:
        sp.check_call(cmd1, stdout=f)

    # call 'denovo.py' through parallelisation wrapper, output to file
    vcf_out_path2 = tempdir.join('out2.vcf').strpath
    cmd2 = [denovo_script, '--parallel',
            '--proband', 'Child',
            '--mother', 'Mother',
            '--father', 'Father',
            datasets['vcf_spw'], vcf_out_path2]

    sp.check_call(cmd2)

    # check that both calls yield the same results
    vr1 = cyvcf2.VCF(vcf_out_path1)
    vr2 = cyvcf2.VCF(vcf_out_path2)

    # same headers
    check_denovo_header(vr1)
    check_denovo_header(vr2)

    # explicitly check that the denovo header lines are the same
    assert get_denovo_program(vr1) == get_denovo_program(vr2)
    assert 'DQ' in get_formats(vr1)

    for variant, variant2 in zip(vr1, vr2):
        assert variants_match(variant, variant2)


def test_cmd_ids_spw_sv_ref(datasets, tempdir):
    """Test cmd interface against SV reference values"""

    cmd1 = [python_exec, denovo_script, '--proband', 'NA12882', '--mother', 'NA12877',
            '--father', 'NA12878', datasets['vcf_spw_sv_chr22']]

    vcf_out_path1 = tempdir.join('.vcf').strpath
    with open(vcf_out_path1, 'w') as f:
        sp.check_call(cmd1, stdout=f)

    vr1 = cyvcf2.VCF(vcf_out_path1)

    check_denovo_header(vr1)

    n_hit = 0
    for variant in vr1:
        assert is_sv(variant)
        v = (variant.CHROM, variant.POS, variant.INFO['SVTYPE'])
        if v in manta_truth_set:
            assert sample_field(variant, 0, 'DQ') >= 0.0
            n_hit += 1

    # check if all truth set variants were called: full recall
    assert n_hit == len(manta_truth_set)


def test_cmd_allosomes_male_proband_spw(datasets, tempdir):
    """Test cmd interface with male-specific variants"""

    cmd1 = [python_exec, denovo_script, '--proband', 'NA12882', '--mother', 'NA12878',
            '--father', 'NA12877', datasets['vcf_spw_chrXY']]

    vcf_out_path1 = tempdir.join('out.vcf').strpath
    with open(vcf_out_path1, 'w') as f:
        sp.check_call(cmd1, stdout=f)

    vr1 = cyvcf2.VCF(vcf_out_path1)

    check_denovo_header(vr1)

    n_denovo = n_no_denovo = n_no_dq = 0
    for variant in vr1:
        # compare against simulation tags
        has_dq_field = 'DQ' in variant.FORMAT
        if not has_dq_field:
            # failing filter in mother or proband
            assert variant.ID == 'noDQ'
            n_no_dq += 1
        else:
            # all other have a DQ score
            if variant.ID == 'denovo':
                assert sample_field(variant, 0, 'DQ') >= 7.0
                n_denovo += 1
            elif variant.ID == 'noDenovo':
                assert sample_field(variant, 0, 'DQ') == 0.0
                n_no_denovo += 1
            elif variant.ID == 'noDQ':
                n_no_dq += 1
            else:
                raise ValueError("No valid annotation for variant")
            # parents always have a missing DQ score
            for i in [1, 2]:
                assert sample_field(variant, i, 'DQ') is None

    assert n_no_dq == 17
    assert n_no_denovo == 7
    assert n_denovo == 5


def test_cmd_pedphase_filter_spw(datasets, tempdir):
    """Test cmd interface with phasing information"""

    cmd1 = [python_exec, denovo_script, '--proband', 'NA12878', '--mother', 'NA12892',
            '--father', 'NA12891', datasets['vcf_spw_pedphase']]

    vcf_out_path1 = tempdir.join('out1.vcf').strpath
    with open(vcf_out_path1, 'w') as f:
        sp.check_call(cmd1, stdout=f)

    cmd2 = cmd1 + ['--ignore-pedigree-phasing']

    vcf_out_path2 = tempdir.join('out2.vcf').strpath
    with open(vcf_out_path2, 'w') as f:
        sp.check_call(cmd2, stdout=f)

    vr1 = cyvcf2.VCF(vcf_out_path1)
    vr2 = cyvcf2.VCF(vcf_out_path2)

    check_denovo_header(vr1)
    check_denovo_header(vr2)

    n_denovo = n_no_denovo = n_no_dq = 0
    for variant, variant2 in zip(vr1, vr2):
        assert variants_match(variant, variant2)
        # compare against simulation tags
        has_dq_field = 'DQ' in variant.FORMAT
        if not has_dq_field:
            # failing filter in mother or proband
            assert variant.ID == 'noDQ'
            n_no_dq += 1
        else:
            # all other have a DQ score
            if variant.ID == 'denovo':
                assert sample_field(variant, 0, 'DQ') >= 0.0
                n_denovo += 1
            elif variant.ID == 'noDenovo':
                assert sample_field(variant, 0, 'DQ') == 0.0
                n_no_denovo += 1
            elif variant.ID == 'noDQ':
                raise ValueError("Variant should not have a DQ score")
            else:
                raise ValueError("No valid annotation for variant")
            # parents always have a missing DQ score
            for i in [1, 2]:
                assert sample_field(variant, i, 'DQ') is None

    assert n_no_dq == 0
    assert n_no_denovo == 7
    assert n_denovo == 4


def test_cmd_proband_sibling_five_samples(datasets, tempdir):
    """Test DQ scores with proband and sibling in two runs [SPW-994]"""

    # VCF samples: 0: Father, 1: Mother, 2: Child, 3: Father2, 4: Mother2

    vcf_proband_path = tempdir.join('proband.vcf').strpath
    vcf_sibling_path = tempdir.join('sibling.vcf').strpath

    cmd_proband = [python_exec, denovo_script,
                   '--proband', 'Child', '--mother', 'Mother', '--father', 'Father',
                   datasets['vcf_spw_5samples']]

    cmd_sibling = [python_exec, denovo_script,
                   '--proband', 'Mother2', '--mother', 'Mother', '--father', 'Father',
                   vcf_proband_path]

    with open(vcf_proband_path, 'w') as f:
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            sp.check_call(cmd_proband, stdout=f)
        assert len(record) == 0  # no warning raised

    with open(vcf_sibling_path, 'w') as f:
        with warnings.catch_warnings(record=True) as record:
            warnings.simplefilter("always")
            sp.check_call(cmd_sibling, stdout=f)
        assert len(record) == 0  # no warning raised

    vf1 = cyvcf2.VCF(vcf_proband_path)
    vf2 = cyvcf2.VCF(vcf_sibling_path)

    check_denovo_header(vf1)
    check_denovo_header(vf2)

    get_dq = lambda variant, sample: sample_field(variant, sample, 'DQ')

    has_dq = lambda variant, sample: has_field(variant, 'DQ')

    # cyvcf2 reports an empty filter set (PASS or '.') as None,
    # mirroring denovo's own can_compute check (FILTER is not None)
    pass_filter = lambda variant, samples: \
                  all([sample_field(variant, i, 'FT')=='PASS' for i in samples]) and \
                  variant.FILTER is None

    sample_indices = lambda vcf, samples: \
                          [vcf.samples.index(sample) for sample in samples]

    for variant1, variant2 in zip(vf1, vf2):
        # compare common most fields are the same
        for field in ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'FILTER', 'FORMAT'):
            assert variant1.__getattribute__(field) == variant2.__getattribute__(field)
        assert info_eq(variant1.INFO, variant2.INFO)
        # compare against simulation tags

        # "proband"
        sample_index = vf1.samples.index('Child')
        assert get_dq(variant1, sample_index) == get_dq(variant2, sample_index)
        trio_indices = sample_indices(vf1, ['Child', 'Mother', 'Father'])
        if pass_filter(variant1, trio_indices):
            assert get_dq(variant1, sample_index) is not None

        # "sibling"
        sample_index = vf1.samples.index('Mother2')
        assert get_dq(variant1, sample_index) is None
        trio_indices = sample_indices(vf2, ['Mother2', 'Mother', 'Father'])
        if pass_filter(variant2, trio_indices):
            assert get_dq(variant2, sample_index) is not None
        else:
            assert get_dq(variant1, sample_index) == get_dq(variant2, sample_index)

        # "unused sample"
        sample_index = vf1.samples.index('Father2')
        assert get_dq(variant1, sample_index) is None
        assert get_dq(variant2, sample_index) is None

        # "mother"
        sample_index = vf1.samples.index('Mother')
        assert get_dq(variant1, sample_index) is None
        assert get_dq(variant2, sample_index) is None

        # "father"
        sample_index = vf1.samples.index('Father')
        assert get_dq(variant1, sample_index) is None
        assert get_dq(variant2, sample_index) is None


def test_cmd_sample_missing(datasets):
    """Test cmd with missing sample"""

    # sample name not present in VCF
    cmd = [python_exec, denovo_script, '--proband', 'MISSING', '--mother',
           'PARENT1', '--father', 'PARENT2', datasets['vcf_in']]

    with pytest.raises(sp.CalledProcessError) as e:
        sp.check_call(cmd, stdout=dev_null, stderr=dev_null)
    assert e.value.returncode > 0

    # incomplete pedigree
    cmd = [python_exec, denovo_script, '--mother', 'PARENT1', '--father',
           'PARENT2', datasets['vcf_in']]

    with pytest.raises(sp.CalledProcessError) as e:
        sp.check_call(cmd, stdout=dev_null, stderr=dev_null)
    assert e.value.returncode > 0


def test_cmd_pedigree_missing(datasets):
    """Test cmd with missing pedigree"""

    cmd = [python_exec, denovo_script, datasets['vcf_in']]

    with pytest.raises(sp.CalledProcessError) as e:
        sp.check_call(cmd, stdout=dev_null, stderr=dev_null)
    assert e.value.returncode > 0


# python interface

def test_python_full(datasets, tempdir):
    """Test python interface against reference values"""

    vcf_out_path = tempdir.join('out1.vcf').strpath
    vcf_out_path2 = tempdir.join('out2.vcf').strpath
    param = {
        'model': 'dng',
        'pass_sample_filters': True,
        'min_read_depth': 1,
        'select_method': 'call'
    }
    pedigree = {
        'proband': 'PROBAND',
        'mother': 'PARENT1',
        'father': 'PARENT2'
    }

    denovo.denovo(datasets['vcf_in'], pedigree, param, vcf_out_path)

    # test if output VCF has the DQ field set
    vr = cyvcf2.VCF(vcf_out_path)

    check_denovo_header(vr)

    n_var1 = 0
    n_dq1 = 0
    for variant in vr:
        assert 'DQ' in variant.FORMAT
        if sample_field(variant, 0, 'DQ') is not None:
            assert sample_field(variant, 0, 'DQ') >= 0.0
            assert sample_field(variant, 0, 'DQ') <= denovo._max_score
            n_dq1 += 1
        assert sample_field(variant, 1, 'DQ') is None
        assert sample_field(variant, 2, 'DQ') is None
        n_var1 += 1

    # test that overwriting the DQ field does not raise a warning
    with warnings.catch_warnings(record=True) as record:
        warnings.simplefilter("always")
        denovo.denovo(vcf_out_path, pedigree, param, vcf_out_path2)
    assert len(record) == 0

    vr2 = cyvcf2.VCF(vcf_out_path2)

    n_var2 = 0
    n_dq2 = 0
    for variant in vr2:
        assert 'DQ' in variant.FORMAT
        if sample_field(variant, 0, 'DQ') is not None:
            assert sample_field(variant, 0, 'DQ') >= 0.0
            assert sample_field(variant, 0, 'DQ') <= denovo._max_score
            n_dq2 += 1
        assert sample_field(variant, 1, 'DQ') is None
        assert sample_field(variant, 2, 'DQ') is None
        n_var2 += 1

    assert n_var1 == n_var2
    assert n_dq1 == n_dq2


def test_python_dummy_model(datasets, tempdir):
    """Test the dummy scoring model via the python interface"""

    vcf_out_path = tempdir.join('out1.vcf').strpath
    param = {
        'model': 'dummy',
        'pass_sample_filters': True,
        'min_read_depth': 1,
        'select_method': 'call'
    }
    pedigree = {
        'proband': 'PROBAND',
        'mother': 'PARENT1',
        'father': 'PARENT2'
    }

    denovo.denovo(datasets['vcf_in'], pedigree, param, vcf_out_path)

    vr = cyvcf2.VCF(vcf_out_path)

    check_denovo_header(vr)

    n_var = 0
    for variant in vr:
        assert 'DQ' in variant.FORMAT
        # dummy model: proband only, constant zero score
        assert sample_field(variant, 0, 'DQ') == 0.0
        assert sample_field(variant, 1, 'DQ') is None
        assert sample_field(variant, 2, 'DQ') is None
        n_var += 1

    assert n_var > 0


def test_python_allosomes_male_proband(datasets, tempdir):
    """Test the python interface for X/Y variants with a male (XY) proband"""

    vcf_out_path = tempdir.join('out1.vcf').strpath
    # the python interface does not merge partial param dicts with the
    # defaults: pass min_gqx explicitly (cli default)
    param = {
        'model': 'dng',
        'pass_sample_filters': True,
        'min_read_depth': 1,
        'select_method': 'call',
        'min_gqx': denovo._min_gqx
    }
    pedigree = {
        'proband': 'NA12882',
        'mother': 'NA12878',
        'father': 'NA12877'
    }

    denovo.denovo(datasets['vcf_spw_chrXY'], pedigree, param, vcf_out_path)

    vr = cyvcf2.VCF(vcf_out_path)

    check_denovo_header(vr)

    n_var = 0
    n_allosomal = 0
    for variant in vr:
        assert 'DQ' in variant.FORMAT
        if variant.CHROM in ('chrX', 'chrY'):
            n_allosomal += 1
            dq = sample_field(variant, 0, 'DQ')
            if dq is not None:
                assert 0.0 <= dq <= denovo._max_score
        assert sample_field(variant, 1, 'DQ') is None
        assert sample_field(variant, 2, 'DQ') is None
        n_var += 1

    assert n_var > 0
    assert n_allosomal > 0


def test_python_sv(datasets, tempdir):
    """Test the python interface for structural variants"""

    vcf_out_path = tempdir.join('out1.vcf').strpath
    param = {
        'model': 'dng',
        'pass_sample_filters': True,
        'min_read_depth': 1,
        'select_method': 'call'
    }
    pedigree = {
        'proband': 'NA12882',
        'mother': 'NA12878',
        'father': 'NA12877'
    }

    denovo.denovo(datasets['vcf_spw_sv_chr22'], pedigree, param, vcf_out_path)

    vr = cyvcf2.VCF(vcf_out_path)

    check_denovo_header(vr)

    # the input header declares DQ as Integer: the output must re-declare it
    # as Float (DQ scores are written with decimals)
    for line in vr.raw_header.splitlines():
        if line.startswith('##FORMAT=') and '<ID=DQ,' in line:
            assert 'Type=Float' in line

    n_var = 0
    n_sv = 0
    for variant in vr:
        assert 'DQ' in variant.FORMAT
        if is_sv(variant):
            n_sv += 1
            dq = sample_field(variant, 0, 'DQ')
            if dq is not None:
                assert 0.0 <= dq <= denovo._max_score
        n_var += 1

    assert n_var > 0
    assert n_sv > 0


def test_python_sample_missing(datasets):
    """Test python interface with missing sample"""

    # sample not in VCF
    pedigree = {
        'proband': 'MISSING',
        'mother': 'PARENT1',
        'father': 'PARENT2'
    }

    with pytest.raises(ValueError):
        denovo.denovo(datasets['vcf_in'], pedigree)

    # sample association missing
    pedigree = {
        'mother': 'PARENT1',
        'father': 'PARENT2'
    }

    with pytest.raises(ValueError):
        denovo.denovo(datasets['vcf_in'], pedigree)


def test_python_read_depth_0(datasets):
    """Test python interface with ill-defined read depth filter"""

    pedigree = {
        'proband': 'PROBAND',
        'mother': 'PARENT1',
        'father': 'PARENT2'
    }
    param = {'model': 'dng', 'min_read_depth': 0}

    with pytest.raises(ValueError):
        denovo.denovo(datasets['vcf_in'], pedigree, param)


def test_python_best_dng_reference(datasets, tempdir):
    """Test python interface with DNG best selection against reference"""

    vcf_out_path1 = tempdir.join('out1.vcf').strpath
    vcf_out_path2 = tempdir.join('out2.vcf').strpath
    param = {
        'model': 'dng',
        'min_read_depth': 1,
        'pass_sample_filters': True,
        'sample_index': [0, 1, 2],
        'select_method': 'best'
    }
    pedigree = {
        'proband': 'PROBAND',
        'mother': 'PARENT1',
        'father': 'PARENT2'
    }

    # unphased
    denovo.denovo(datasets['vcf_in'], pedigree, param, vcf_out_path1)
    vr1 = cyvcf2.VCF(vcf_out_path1)

    # phased
    denovo.denovo(datasets['vcf_in_phased'], pedigree, param, vcf_out_path2)
    vr2 = cyvcf2.VCF(vcf_out_path2)

    # reference
    ref = gzip.open(datasets['dng_ref'])

    check_denovo_header(vr1)
    check_denovo_header(vr2)

    for v1, v2, r in zip(vr1, vr2, ref):
        rf = r.split(b' ')
        # same position
        assert v1.POS == int(rf[6])
        assert v2.POS == int(rf[6])
        # compare denovo scores with reference values
        dng = round(denovo.prob2score(1.0 - float(rf[24])), 1)
        if 'DQ' in v1.FORMAT and sample_field(v1, 0, 'DQ') is not None:
            assert approx(sample_field(v1, 0, 'DQ'), dng)
        if 'DQ' in v2.FORMAT and sample_field(v2, 0, 'DQ') is not None:
            assert approx(sample_field(v2, 0, 'DQ'), dng)


def test_python_call_dng_reference(datasets, tempdir):
    """Test python interface with DNG call selection against reference"""

    vcf_out_path = tempdir.join('out.vcf').strpath
    param = {
        'model': 'dng',
        'min_read_depth': 1,
        'pass_sample_filters': True,
        'sample_index': [0, 1, 2],
        'select_method': 'call'
    }
    pedigree = {
        'proband': 'PROBAND',
        'mother': 'PARENT1',
        'father': 'PARENT2'
    }

    denovo.denovo(datasets['vcf_in'], pedigree, param, vcf_out_path)

    # test if output VCF has the DQ field set
    ref = gzip.open(datasets['dng_ref'])
    vr = cyvcf2.VCF(vcf_out_path)

    check_denovo_header(vr)

    n_hits = 0
    for variant, r in zip(vr, ref):
        rf = r.split(b' ')
        # same position
        assert variant.POS == int(rf[6])
        # compare denovo scores with reference values
        dng = round(denovo.prob2score(1.0 - float(rf[24])), 1)
        if 'DQ' in variant.FORMAT:
            dq = sample_field(variant, 0, 'DQ')
            # py2 semantics: a missing (None) score is trivially <= dng
            if dq is not None:
                assert dq <= float(dng) + approx(0.0, 1.0)
            n_hits += 1

    assert n_hits > 0


# SPW/strelka multi-sample input

def test_python_spw_ref(datasets, tempdir):
    """Test python interface for SPW data against reference values"""

    vcf_out_path = tempdir.join('out.vcf').strpath
    param = {'model': 'dng',
             'min_read_depth': 10,
             'pass_sample_filters': True,
             'select_method': 'call'}
    pedigree = {'proband': 'Child', 'mother': 'Mother', 'father': 'Father'}

    denovo.denovo(datasets['vcf_spw'], pedigree, param, vcf_out_path)

    vr = cyvcf2.VCF(vcf_out_path)
    child_idx = vr.samples.index(pedigree['proband'])

    check_denovo_header(vr)

    n_hit = 0
    n_called = {'indel': 0, 'snp': 0}
    for variant in vr:
        v = (variant.CHROM, variant.POS, variant.REF, variant.ALT)
        if 'DQ' in variant.FORMAT:
            dq = sample_field(variant, child_idx, 'DQ')
            if dq is not None and dq >= 15.0:
                n_called[var_type(variant)] += 1
            for sample_index in range(len(variant.genotypes)):
                if sample_index != child_idx:
                    assert sample_field(variant, sample_index, 'DQ') is None
                else:
                    assert sample_field(variant, sample_index, 'DQ') is None or sample_field(variant, sample_index, 'DQ') >= 0.0
            # check number of decimals
            dq = sample_field(variant, child_idx, 'DQ')
            if dq is not None:
                assert approx(dq, round(dq, denovo._default_dq_decimals))
                assert dq <= denovo._max_score
        if v in strelka_truth_set:
            assert sample_field(variant, child_idx, 'DQ') >= 20.0
            n_hit += 1

    assert n_called['snp'] >= 4
    assert n_called['indel'] >= 2
    # check if all truth set variants were called: full recall
    assert n_hit == len(strelka_truth_set)


def test_python_spw_decimals_ref(datasets, tempdir):
    """Test python interface for SPW data with custom number of decimals"""

    vcf_out_path = tempdir.join('out.vcf').strpath
    decimals = 5
    param = {
        'model': 'dng',
        'min_read_depth': 10,
        'pass_sample_filters': True,
        'select_method': 'call',
        'score_decimals': decimals}
    pedigree = {'proband': 'Child', 'mother': 'Mother', 'father': 'Father'}

    denovo.denovo(datasets['vcf_spw'], pedigree, param, vcf_out_path)

    vr = cyvcf2.VCF(vcf_out_path)
    child_idx = vr.samples.index(pedigree['proband'])

    check_denovo_header(vr)

    n_hit = 0
    n_called = {'indel': 0, 'snp': 0}
    for variant in vr:
        assert has_field(variant, 'DQ')
        v = (variant.CHROM, variant.POS, variant.REF, variant.ALT)
        if 'DQ' in variant.FORMAT:
            dq = sample_field(variant, child_idx, 'DQ')
            if dq is not None and dq >= 15.0:
                n_called[var_type(variant)] += 1
            for sample_index in range(len(variant.genotypes)):
                if sample_index != child_idx:
                    assert sample_field(variant, sample_index, 'DQ') is None
                else:
                    sample_field(variant, sample_index, 'DQ') is None or sample_field(variant, sample_index, 'DQ') >= 0.0
            # check number of decimals
            dq = sample_field(variant, child_idx, 'DQ')
            if dq is not None:
                assert approx(dq, round(dq, decimals))
                assert dq <= denovo._max_score
        if v in strelka_truth_set:
            assert sample_field(variant, child_idx, 'DQ') >= 20.0
            n_hit += 1

    assert n_called['snp'] >= 4
    assert n_called['indel'] >= 2
    # check if all truth set variants were called: full recall
    assert n_hit == len(strelka_truth_set)


def test_python_spw_ref_five_samples(datasets, tempdir):
    """Test python interface for 5 sample SPW data against reference values"""

    vcf_out_path = tempdir.join('out.vcf').strpath
    param = {'model': 'dng',
             'min_read_depth': 10,
             'pass_sample_filters': True,
             'select_method': 'call'}
    pedigree = {'proband': 'Child', 'mother': 'Mother', 'father': 'Father'}

    denovo.denovo(datasets['vcf_spw_5samples'], pedigree, param, vcf_out_path)

    vr = cyvcf2.VCF(vcf_out_path)
    child_idx = vr.samples.index(pedigree['proband'])

    check_denovo_header(vr)

    n_hit = 0
    n_called = {'indel': 0, 'snp': 0}
    for variant in vr:
        v = (variant.CHROM, variant.POS, variant.REF, variant.ALT)
        if v in strelka_truth_set:
            assert sample_field(variant, child_idx, 'DQ') >= 20.0
            n_hit += 1
        if 'DQ' in variant.FORMAT:
            for sample_index in range(len(variant.genotypes)):
                if sample_index != child_idx:
                    assert sample_field(variant, sample_index, 'DQ') is None
                elif (dq := sample_field(variant, sample_index, 'DQ')) is not None and dq >= 15.0:
                    n_called[var_type(variant)] += 1

    assert n_called['snp'] >= 4
    assert n_called['indel'] >= 2
    # check if all truth set variants were called: full recall
    assert n_hit == len(strelka_truth_set)


def test_python_spw_old_prior_cases(datasets, tempdir):
    """Test cases in which the old prior prevented us from calling denovo events"""

    vcf_out_path = tempdir.join('out.vcf').strpath
    param = {'model': 'dng',
             'min_read_depth': 10,
             'pass_sample_filters': True,
             'select_method': 'call'}
    pedigree = {'proband': 'NA12877', 'mother': 'NA12878', 'father': 'NA12882'}

    denovo.denovo(datasets['vcf_spw1005'], pedigree, param, vcf_out_path)

    vr = cyvcf2.VCF(vcf_out_path)
    assert 'DQ' in get_formats(vr)

    n_denovo = n_no_denovo = n_no_dq = 0
    for variant in vr:
        # compare against simulation tags
        assert 'DQ' in variant.FORMAT
        # all other have a DQ score
        if variant.ID == 'denovo':
            assert sample_field(variant, 0, 'DQ') >= 7.0
            n_denovo += 1
        elif variant.ID == 'noDenovo':
            assert sample_field(variant, 0, 'DQ') < 2.0
            n_no_denovo += 1
        else:
            raise ValueError("No valid annotation for variant")
        # parents always have a missing DQ score
        for i in [1, 2]:
            assert sample_field(variant, i, 'DQ') is None

    assert n_no_dq == 0
    assert n_no_denovo == 4
    assert n_denovo == 2


# prior import, reordering and computation

def test_snv_prior_import_order():
    """Test SNV prior import and reordering"""

    for prior_type, prior_path in denovo._prior_paths.items():

        prior_path = prior_path['snv']
        assert os.path.exists(prior_path)

        # import
        prior = denovo.read_prior(prior_path, 'snv')
        assert len(prior) == 11  # number of columns
        assert len(prior['gt']) == 1000  # number of rows
        assert len(set(prior['gt'])) == 1000  # only unique genotypes

        # reordering
        prior = denovo.reorder_prior(prior, denovo._full_gt_idx_snv)
        # check order of genotypes
        gt = prior['gt']
        assert gt[0] == 'AA/AA/AA'
        assert gt[15] == 'CG/AC/AA'
        assert gt[998] == 'GT/TT/TT'
        assert gt[999] == 'TT/TT/TT'


def test_indel_length_prior_computation():
    """Test computation of length-dependent indel prior"""

    max_len = 5
    p = {'slope': -0.5, 'intercept': -1.5}
    mut_rate = denovo.precompute_indel_prior(p, max_len)
    assert len(mut_rate) == max_len+1
    assert np.all(np.diff(mut_rate) < 0.0)  # decaying values

    for prior_type, prior_path in denovo._prior_paths.items():

        prior_path = prior_path['indel']
        assert os.path.exists(prior_path)

        # import
        prior = denovo.read_prior(prior_path, 'indel')
        assert len(prior) == 10  # number of columns
        assert len(prior['gt']) == 27  # number of rows
        assert len(set(prior['gt'])) == 27  # only unique genotypes

        # reordering
        prior = denovo.reorder_prior(prior, denovo._full_gt_idx_indel)
        # check order of genotypes
        gt = prior['gt']
        assert gt[0] == 'RR/RR/RR'
        assert gt[1] == 'RV/RR/RR'
        assert gt[25] == 'RV/VV/VV'
        assert gt[26] == 'VV/VV/VV'

        for indel_type in ('ins', 'del'):
            p = prior[indel_type]
            assert p.shape == (27, denovo._max_indel_len+1)
            for i in range(27):
                assert np.all(np.diff(p[i, :]) <= 0.0)  # decaying values
                # TODO check for strictly monotonic decreasing values: i=5


def test_mutation_at_reference_lowlevel(datasets):
    """Test python interface against DNG reference values"""

    vcf_path_full = datasets['mutation_AT']
    vcf_path_min = datasets['mutation_AT_minimal']

    for test_vcf_path in [vcf_path_full, vcf_path_min]:

        snv_prior = denovo.read_prior(denovo._prior_paths['auto']['snv'], 'snv')
        snv_prior = denovo.reorder_prior(snv_prior, denovo._full_gt_idx_snv)
        prior = {'auto': {'snv': snv_prior, 'indel': None}}

        vf = cyvcf2.VCF(test_vcf_path)
        v = next(vf)

        DQ_score_ref = denovo._max_score

        param = {
            'min_read_depth': 1,
            'pass_sample_filters': False,
            'is_select_call': True,
            'filter_pedphase': False
        }

        # correct trio
        param['sample_index'] = [0, 1, 2]
        dq1 = denovo.calculate_dng_DQ(v, prior, param)
        assert round(dq1, 3) == DQ_score_ref

        # switching mother and father should have no effect
        param['sample_index'] = [0, 2, 1]
        dq2 = denovo.calculate_dng_DQ(v, prior, param)
        assert round(dq1, 10) == round(dq2, 10)

        # switching the child with any parent should have an effect
        param['sample_index'] = [1, 0, 2]
        dq2 = denovo.calculate_dng_DQ(v, prior, param)
        assert round(dq1, 10) != round(dq2, 10)

        # switching the child with any parent should have an effect
        param['sample_index'] = [2, 1, 0]
        dq2 = denovo.calculate_dng_DQ(v, prior, param)
        assert round(dq1, 10) != round(dq2, 10)

        # 'best' select method
        param['sample_index'] = [0, 1, 2]
        param['select_method'] = 'best'
        dq3 = denovo.calculate_dng_DQ(v, prior, param)
        assert round(dq3, 3) == DQ_score_ref


def test_vcf_missing_essential_fields(datasets):
    """Test exception handling of VCF without essential fields"""

    vcf_path2 = datasets['missing_PL']

    pedigree = {
        'proband': 'child',
        'father': 'father',
        'mother': 'mother'
    }

    for vcf_path in [vcf_path2]:
        with pytest.raises(ValueError):
            denovo.denovo(vcf_path, pedigree, denovo._default_params)


# score conversion

def test_score_conversion():
    """Test conversion between score and probabilities"""

    atol = 1e-5
    ms = denovo._max_score
    prob = [0.0,  1e-100, 1e-10, 0.01, 0.1,  0.794328, 0.5,    1.0]
    score = [ms,  ms,     100.0, 20.0, 10.0, 1.0,      3.0103, 0.0]

    assert len(prob) == len(score)

    for p, q in zip(prob, score):
        assert_allclose(denovo.prob2score(p), q, atol=atol)
        assert_allclose(denovo.score2prob(q), p, atol=atol)

    assert_allclose(denovo.prob2score_vec(prob), score, atol=atol)
    assert_allclose(denovo.score2prob_vec(score), prob, atol=atol)


def test_allele_snv_indices():
    """Test allele indices calculation for SNVs"""

    bases = ['A', 'G', 'C']
    combs = [('A', 'A'), ('A', 'G'), ('G', 'G'),
             ('A', 'C'), ('G', 'C'), ('C', 'C')]

    idx = denovo.allele_indices(bases, denovo._gts_idx_snv)
    assert idx == [0, 2, 7, 1, 5, 4]
    assert idx == [denovo._gts_idx_snv[x] for x in combs]


def test_allele_indel_indices():
    """Test allele indices calculation for indels"""

    bases = ['R', 'V', 'V']
    idx = denovo.allele_indices(bases, denovo._gts_idx_indel)
    assert idx == [0, 1, 2, 1, 2, 2]


def test_PL_indices():
    """Test alleles to indices calculation for PL fields"""

    # generate lookup tables
    pl_idx_snv = denovo.build_PL_indices_lookup(
        denovo._alleles_snv, denovo._gts_idx_snv)
    pl_idx_indel = denovo.build_PL_indices_lookup(
        denovo._alleles_indel, denovo._gts_idx_indel, 'R')

    assert len(pl_idx_snv) == 60

    assert len(pl_idx_indel) == 3

    assert pl_idx_snv[('A', 'C')] == [0, 1, 4]
    assert pl_idx_snv[('C', 'A')] == [4, 1, 0]
    assert pl_idx_snv[('A', 'C', 'G', 'T')] == [0, 1, 4, 2, 5, 7, 3, 6, 8, 9]
    assert pl_idx_snv[('G', 'A', 'T', 'C')] == [7, 2, 0, 8, 3, 9, 5, 1, 6, 4]

    assert pl_idx_indel[('R', 'V')] == [0, 1, 2]
    assert pl_idx_indel[('R', 'V', 'V')] == [0, 1, 2, 1, 2, 2]
    assert pl_idx_indel[('R', 'V', 'V', 'V')] == [0, 1, 2, 1, 2, 2, 1, 2, 2, 2]

    # lookup
    assert denovo.get_PL_indices_snv(('A', 'C')) == [0, 1, 4]
    assert denovo.get_PL_indices_indel(('R', 'V')) == [0, 1, 2]
    # freshly computed
    assert len(denovo.get_PL_indices_snv(('A', 'C', 'A', 'C', 'A', 'C'))) == 21
    assert len(denovo.get_PL_indices_indel(('V', 'V', 'V', 'V', 'V', 'V'))) == 21


def test_utils_exceptions(datasets):
    """Test exception handling of low-level functions"""

    # non-existing data file
    with pytest.raises(IOError):
        denovo.get_data_path('nofile', 'nodir')

    with pytest.raises(IOError):
        get_test_data_path('nofile', 'nodir')

    # unknown variant type for prior
    with pytest.raises(ValueError):
        denovo.read_prior('some_path', 'SV')

    # input VCF missing
    pedigree = {'proband': 'Child', 'mother': 'Mother', 'father': 'Father'}
    with pytest.raises(IOError):
        denovo.denovo('no_file', pedigree)

    # select method not valid
    param = denovo._default_params
    param['select_method'] = 'invalid'
    param['model'] = denovo._default_model_name
    with pytest.raises(ValueError):
        denovo.denovo(datasets['vcf_spw'], pedigree, param)


def test_denovo_invalid_vcf_input(datasets):
    """Test exception handling for invalid VCF file"""

    # input some non-VCF file
    pedigree = {
        'proband': 'PROBAND',
        'mother': 'PARENT2',
        'father': 'PARENT1'
    }

    with pytest.raises(IOError):
        denovo.denovo(datasets['invalid_vcf'], pedigree)


def test_denovo_catch_exceptions(datasets, tempdir):
    """Test catching of potential exceptions"""

    vcf_out_path = tempdir.join('out.vcf').strpath

    cmd = [python_exec, denovo_script, '--proband', 'CephChild-12882',
           '--mother', 'CephMother-12878', '--father', 'CephFather-12877',
           datasets['vcf_spw542']]

    # No exceptions raised at user-level
    with open(vcf_out_path, 'w') as f:
        sp.check_call(cmd, stdout=f)


def test_pedigree_parsing_from_vcf_header(datasets):
    """Test extraction of pedigree from SPW VCF"""

    # VCF file without pedigree information
    vf = cyvcf2.VCF(datasets['vcf_spw'])
    assert denovo.parse_pedigree_from_vcf_header(vf) is None

    # SPW VCF file with pedigree information
    vf = cyvcf2.VCF(datasets['vcf_spw_chr22'])
    pedigree = denovo.parse_pedigree_from_vcf_header(vf)
    pedigree_check = denovo.check_pedigree_input(pedigree, vf.samples)
    assert pedigree == pedigree_check


def test_karyotypes_parsing_from_vcf_header(datasets):
    """Test extraction of predicted karyotypes from SPW VCF"""

    # VCF file without pedigree information
    vf = cyvcf2.VCF(datasets['vcf_spw'])
    assert denovo.parse_karyotypes_from_vcf_header(vf) is None

    # SPW VCF file with pedigree information
    ref_karyo = {
        'CephFather-12877': 'XY',
        'CephMother-12878': 'XX',
        'CephChild-12882': 'XY'
    }
    vf = cyvcf2.VCF(datasets['vcf_spw_chr22'])
    karyotypes = denovo.parse_karyotypes_from_vcf_header(vf)
    assert ref_karyo == karyotypes


def test_cyvcf2_variant_examples(datasets):
    """Test basic cyvcf2 accessors and variant types"""

    example_reference = [
        ('chr20', 60000000, 'A', ['GCA'], 'indel'),
        ('chr20', 60000001, 'ACG', ['T'], 'indel',),
        ('chr20', 60000002, 'A', ['G', 'GCA'], 'indel'),
        ('chr20', 60000003, 'A', ['G'], 'snp'),
        ('chr20', 60000004, 'A', ['G', 'T'], 'snp')
    ]

    vf = cyvcf2.VCF(datasets['variant_examples'])

    for v, r in zip(vf, example_reference):
        assert v.CHROM == r[0]
        assert v.POS == r[1]
        assert v.REF == r[2]
        assert v.ALT == r[3]
        assert var_type(v) == r[4]
        if r[4] == 'snp':
            assert is_snp(v)
        elif r[4] == 'indel':
            assert is_indel(v)


def test_cyvcf2_exceptions(datasets):
    """Test exception handling of cyvcf2 for invalid VCF inputs"""

    with pytest.raises(Exception):
        cyvcf2.VCF("no_vcf_file")

    with pytest.raises(IOError):
        cyvcf2.VCF(datasets['invalid_vcf'])


def test_indel_length(datasets):
    """Test computation of indel length"""

    vcf_example_path = datasets['variant_examples']

    example_reference = [
        ('chr20', 60000000, 'A', ['GCA'], 'indel', 2),
        ('chr20', 60000001, 'ACG', ['T'], 'indel', 2),
        ('chr20', 60000002, 'A', ['G', 'GCA'], 'indel', 0),
        ('chr20', 60000003, 'A', ['G'], 'snp', 0),
        ('chr20', 60000004, 'A', ['G', 'T'], 'snp', 0)
    ]

    vf = cyvcf2.VCF(vcf_example_path)

    for v, r in zip(vf, example_reference):
        assert v.CHROM == r[0]
        assert v.POS == r[1]
        assert v.REF == r[2]
        assert v.ALT == r[3]
        assert var_type(v) == r[4]
        assert denovo.indel_length(v) == r[5]


def test_PL_field_formatting(datasets):
    """Test reformatting of PL field in case of different ploidy across samples"""

    vf = cyvcf2.VCF(datasets['vcf_spw_chrXY'])

    v = next(vf)
    pl_old = v.format('PL')
    pl_new = denovo.rectify_pl_field(pl_old)

    for sample in [0, 1, 2]:
        assert pl_old[sample][0] == pl_new[sample][0]
        assert pl_new[sample].min() >= 0

    for sample in [0, 1, 2]:
        assert pl_old[sample][0] == pl_new[sample][0]
        assert pl_new[sample].min() >= 0
        if pl_old[sample][2] < 0:  # sample 0 and 2
            assert pl_new[sample[2]] == pl_old[sample][1]
        else:  # sample unchanged
            assert all(pl_new[sample] == pl_new[sample])


def test_is_field_in_vcf_header(datasets):
    """Test if VCF fields are properly checked"""

    vf = cyvcf2.VCF(datasets['vcf_spw_pedphase'])

    sections = {'GENERIC', 'CONTIG', 'INFO', 'FORMAT', 'STR'}

    references = (
        ('ME', 'FORMAT'),
        ('chr1', 'CONTIG'),
        ('PL', 'FORMAT'),
        ('LowGQX', 'FILTER')
    )

    for (field, section) in references:
        assert denovo.is_field_in_vcf_header(vf, field, section)
        for other_section in sections.difference([section]):
            assert not denovo.is_field_in_vcf_header(vf, field, other_section)

    for section in sections:
        assert not denovo.is_field_in_vcf_header(vf, 'not_existing', section)

    with pytest.raises(ValueError):
        denovo.is_field_in_vcf_header(vf, 'PL', 'non_existing_section')


def test_is_denovo_candidate_pedphase_mendelconflict(datasets, tempdir):
    """Test if candidate selection works with new pedphase ME format field"""

    vf = cyvcf2.VCF(datasets['vcf_spw_pedphase'])

    param = {
        'filter_pedphase': True,
        'sample_index': [0, 1, 2]  # child is the proband
    }
    param_swapped = {
        'filter_pedphase': True,
        'sample_index': [1, 2, 0]  # a parent is the proband
    }

    assert denovo._mendel_conflict_name in vf

    for variant in vf:
        assert denovo._mendel_conflict_name in variant.FORMAT
        me = variant.format('ME')
        proband_has_mendel_conflict = me[0][0] == 1
        is_candidate = denovo.is_denovo_candidate(variant, param)

        if proband_has_mendel_conflict:
            assert is_candidate
            assert variant.ID == 'denovo'
        else:
            assert not is_candidate
            assert variant.ID == 'noDenovo'

        # check possible values for all samples
        # proband: yes or no
        assert me[0][0] in [0, 1]
        # parents: always missing
        assert me[1][0] < 0
        assert me[2][0] < 0
        # sanity check: an invalid ME field would be treated as a candidate
        is_candidate_swapped = denovo.is_denovo_candidate(variant, param_swapped)
        assert is_candidate_swapped

    # VCF w/o ME field
    vf2 = cyvcf2.VCF(datasets['vcf_spw'])

    assert denovo._mendel_conflict_name not in vf2

    for variant in vf2:
        is_candidate = denovo.is_denovo_candidate(variant, param)
        assert is_candidate

        assert denovo._mendel_conflict_name not in variant.FORMAT
        with pytest.raises(KeyError):
            variant.format('ME')


def test_import_par_bed(datasets):
    """Test the import of the PAR from the annotation BED file"""

    par_path = datasets['par_bed']

    par = denovo.import_bed_regions(par_path)

    # check properties of data structure
    for contig, coords in par.items():
        assert contig == 'X'
        assert isinstance(coords, tuple)
        for pos in coords:
            assert pos[0] < pos[1]
            assert len(pos) == 2

    # check properties of the PARs
    assert len(par) == 1  # one contig
    assert len(par['X']) == 2  # two entries

    # check exception handling
    assert denovo.import_bed_regions('not_a_file') is None
    assert denovo.import_bed_regions(datasets['par_vcf']) is None


def test_is_within_par(datasets):
    """Test the import of the PAR from the annotation BED file"""

    par = denovo.import_bed_regions(datasets['par_bed'])
    vf = cyvcf2.VCF(datasets['par_vcf'])

    for variant in vf:
        is_in_par = variant.ID == 'insidePAR'
        assert denovo.is_in_regions(variant, par) == is_in_par

        # in case also the contig is supplied
        chrom_stripped = variant.CHROM.lstrip('chr')
        assert denovo.is_in_regions(variant, par, chrom_stripped) == denovo.is_in_regions(variant, par)

        assert not denovo.is_in_regions(variant, par, 'not_the_right_contig')

        # in case the regions are missing
        assert not denovo.is_in_regions(variant, None)
        assert not denovo.is_in_regions(variant, None, chrom_stripped)
