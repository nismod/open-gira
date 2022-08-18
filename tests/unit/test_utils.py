import numpy as np
import pandas as pd
import pytest

from workflow.scripts.transport.utils import str_to_bool


class TestStrToBool:
    def test_str_mapping(self):
        truth_table = (
            ('N', False),
            ('n', False),
            ('NO', False),
            ('no', False),
            ('FALSE', False),
            ('false', False),
            ('F', False),
            ('f', False),
            ('', False),
            # N.B. null values are transformed to empty strings
            # then mapped to false
            (None, False),
            (np.nan, False),
            ('any other string', True),
        )
        strings, bools = zip(*truth_table)
        expected = pd.Series(bools)
        assert expected.equals(str_to_bool(pd.Series(strings)))

    def test_incorrect_type(self):
        for data in [True, 1, 1.414]:
            with pytest.raises(ValueError):
                str_to_bool(pd.Series([data]))
