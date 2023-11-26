# import pytest
# import os
#
# @pytest.fixture
# def ec_matrix_file():
#     "creates an ec_file with 10 entries"
#     import tempfile
#     fname = tempfile.mktemp()
#     with open(fname, 'w') as fh:
#         for i in range(10):
#             fh.write(f'{i} 1,2,3,4\n')
#     yield fname
#     os.remove(fname)
#
#
# @pytest.fixture
# def transcript_file():
#     "creates an transcript_file with 10 entries"
#     import tempfile
#     fname = tempfile.mktemp()
#     with open(fname, 'w') as fh:
#         for i in range(10):
#             fh.write(f'ENST00000000000{i}.1\n')
#     yield fname
#     os.remove(fname)
