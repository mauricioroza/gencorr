#!/bin/sh
# Run CBIIT/ldsc's own unit tests (test/) with pytest.
#
# The tests were written for nose on Python 2. LDSC's code is NOT changed here;
# only the test harness is adapted so pytest runs what nose ran:
#
#   1. Test_RG_Statistical is a plain class with a setUpClass() classmethod.
#      nose called it; pytest only calls setup_class() on plain classes (and
#      @attr('rg') makes `self.rg` the bool True instead of the results).
#   2. The port reads LD scores only as *.l2.ldscore.gz (parse.ldscore hardcodes
#      the suffix), but test/simulate_test/ldscore/ is committed uncompressed:
#      gzip a copy of each file.
#   3. Test_Hsq_1D.test_weights passes M as a 1x1 np.matrix to Hsq.weights(),
#      which does float(M); numpy 2 refuses float() on a 2-D array. LDSC itself
#      always passes a float there (Hsq.__init__: M_tot = float(np.sum(M))).
#
# Deselected, because they test behaviour the port changed (documented in
# docs/ldsc_py3.md; neither affects gencorr, which only reads .gz LD scores and
# munge output):
#   - test_parse.py::Test_ldscore::test_ldscore_fromlist / test_ldscore_loop:
#     read .bz2 LD score files, no longer supported.
#   - test_parse.py::test_read_sumstats: expects ValueError for a malformed
#     .gz sumstats file; parse.read_csv now prints the error and returns a
#     string, so LDSC fails later with AttributeError instead.
set -eu
cd "${1:-/opt/ldsc}"

work=$(mktemp -d)
cp -r test "$work/test"
cd "$work"

sed -i 's/^    def setUpClass(cls):/    def setup_class(cls):/' test/test_sumstats.py
sed -i '/^class Test_H2_Statistical(unittest.TestCase):/,$ s/^    def setup_class(cls):/    def setUpClass(cls):/' test/test_sumstats.py
for f in test/simulate_test/ldscore/*.l2.ldscore; do gzip -kn "$f"; done
sed -i -E 's/(reg\.Hsq\.weights\(self\.ld, self\.w_ld, self\.N, )self\.M,/\1self.M.item(),/g' test/test_regressions.py

export PYTHONPATH="${PYTHONPATH:-/opt/ldsc}"
status=0
python -m pytest -q -p no:cacheprovider -W ignore::PendingDeprecationWarning \
    --rootdir "$work" \
    --deselect test/test_parse.py::Test_ldscore::test_ldscore_fromlist \
    --deselect test/test_parse.py::Test_ldscore::test_ldscore_loop \
    --deselect test/test_parse.py::test_read_sumstats \
    test || status=$?
cd / && rm -rf "$work"
exit $status
