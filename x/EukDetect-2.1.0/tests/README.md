# EukDetect Tests

## Installation

### Install EukDetect in development mode:
```bash
# Create conda environment from environment.yml
conda env create -f eukdetect/envs/eukdetect2_environment.yml
conda activate eukdetect

# Install EukDetect in editable mode
pip install -e .
```

The conda environment includes pytest and pytest-cov, so no additional test dependencies are needed.

## Running Tests

### Run all tests:
```bash
# From project root
pytest tests/ -v

# With coverage report
pytest tests/ -v --cov=eukdetect --cov-report=html

# Run specific test file
pytest tests/test_eukdetect.py -v

# Run specific test class
pytest tests/test_eukdetect.py::TestValidateCores -v

# Run specific test
pytest tests/test_eukdetect.py::TestValidateCores::test_valid_cores -v
```

## Test Coverage

Current test coverage:
- `validate_cores()` - ✅ Complete (4 tests)
- `validate_sample_name()` - ✅ Complete (7 tests)
- `_parse_samples()` - ✅ Complete (9 tests)
- Config building - ✅ Complete (3 tests)
- File validation - ✅ Complete (6 tests)
- Snakemake executor - ✅ Complete (3 tests)

**Total: 32 comprehensive tests**

## Test Fixtures

The test suite includes fixtures that create dummy files:
- `test_fastq_files` - Creates `test.fastq.gz`, `test_1.fastq.gz`, `test_2.fastq.gz`
- `test_database` - Creates minimal database structure with all required files

These fixtures use pytest's `tmp_path` for automatic cleanup.

## Writing New Tests

### Test Structure
```python
class TestFeatureName:
    """Test description"""
    
    def test_specific_behavior(self):
        """Test that X does Y"""
        # Arrange
        input_data = "test"
        
        # Act
        result = function_to_test(input_data)
        
        # Assert
        assert result == expected_output
```

### Using Fixtures
```python
def test_with_fixture(test_fastq_files):
    """Pytest provides test_fastq_files fixture automatically"""
    se_file = test_fastq_files['se']
    pe1_file = test_fastq_files['pe1']
    pe2_file = test_fastq_files['pe2']
    # Test with these files
```

### Testing Exceptions
```python
def test_invalid_input_raises_error():
    """Test that function raises appropriate error"""
    with pytest.raises(ValueError, match="error message"):
        function_that_should_fail("bad input")
```

## CI Integration

Add to GitHub Actions workflow:
```yaml
name: Tests
on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: conda-incubator/setup-miniconda@v2
        with:
          environment-file: eukdetect/envs/eukdetect2_environment.yml
          activate-environment: eukdetect
      - name: Install package
        shell: bash -l {0}
        run: pip install -e .
      - name: Run tests
        shell: bash -l {0}
        run: pytest tests/ -v --cov=eukdetect --cov-report=xml
      - name: Upload coverage
        uses: codecov/codecov-action@v3
```

## Test Data Files

Tests use hardcoded filenames as specified:
- **Single-end:** `test.fastq.gz`
- **Paired-end:** `test_1.fastq.gz` and `test_2.fastq.gz`

These are automatically created by the `test_fastq_files` fixture.
