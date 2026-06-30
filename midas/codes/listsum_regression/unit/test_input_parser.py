"""
Unit tests for the input parser module.

Tests the parsing and validation of YAML input files,
including error handling and edge cases.

Author: MIDAS Development Team
"""

import pytest
import yaml
from pathlib import Path
from unittest.mock import patch, mock_open
import tempfile

from midas.input_parser import Input_Parser, yaml_line_reader, validate_input


class TestYamlLineReader:
    """Test the yaml_line_reader function."""
    
    def test_yaml_line_reader_with_existing_key(self):
        """Test reading an existing key from YAML data."""
        data = {'test_key': 'test_value'}
        result = yaml_line_reader(data, 'test_key', 'default_value')
        assert result == 'test_value'
    
    def test_yaml_line_reader_with_missing_key(self):
        """Test reading a missing key returns default value."""
        data = {'other_key': 'other_value'}
        result = yaml_line_reader(data, 'missing_key', 'default_value')
        assert result == 'default_value'
    
    def test_yaml_line_reader_with_empty_data(self):
        """Test reading from empty data returns default."""
        data = {}
        result = yaml_line_reader(data, 'any_key', 'default_value')
        assert result == 'default_value'


class TestValidateInput:
    """Test the validate_input function."""
    
    def test_validate_debug_mode_true(self):
        """Test debug mode validation with True value."""
        result = validate_input('debug_mode', True)
        assert result is True
    
    def test_validate_debug_mode_false(self):
        """Test debug mode validation with False value."""
        result = validate_input('debug_mode', False)
        assert result is False
    
    def test_validate_debug_mode_string(self):
        """Test debug mode validation with string value."""
        result = validate_input('debug_mode', 'true')
        assert result is True
    
    def test_validate_results_directory_name(self):
        """Test results directory name validation."""
        result = validate_input('results_directory_name', 'test dir')
        assert result == Path('test_dir')
    
    def test_validate_set_seed_with_value(self):
        """Test set_seed validation with a value."""
        result = validate_input('set_seed', 12345)
        assert result == 12345
    
    def test_validate_set_seed_none(self):
        """Test set_seed validation with None."""
        result = validate_input('set_seed', None)
        assert result is None
    
    def test_validate_unknown_key(self):
        """Test validation of unknown key returns original value."""
        result = validate_input('unknown_key', 'test_value')
        assert result == 'test_value'


class TestInputParser:
    """Test the Input_Parser class."""
    
    def test_input_parser_initialization(self, sample_yaml_file):
        """Test Input_Parser initialization with valid YAML file."""
        parser = Input_Parser(cpus=2, input_file=str(sample_yaml_file))
        assert parser is not None
        assert parser.cpus == 2
    
    def test_input_parser_with_invalid_file(self, temp_dir):
        """Test Input_Parser with non-existent file."""
        invalid_file = temp_dir / "nonexistent.yaml"
        with pytest.raises(FileNotFoundError):
            Input_Parser(cpus=1, input_file=str(invalid_file))
    
    def test_input_parser_with_invalid_yaml(self, temp_dir):
        """Test Input_Parser with invalid YAML content."""
        invalid_yaml = temp_dir / "invalid.yaml"
        with open(invalid_yaml, 'w') as f:
            f.write("invalid: yaml: content: [")
        
        with pytest.raises(yaml.YAMLError):
            Input_Parser(cpus=1, input_file=str(invalid_yaml))
    
    @patch('midas.input_parser.validate_input')
    def test_input_parser_calls_validation(self, mock_validate, sample_yaml_file):
        """Test that Input_Parser calls validate_input for parsed values."""
        mock_validate.return_value = 'validated_value'
        parser = Input_Parser(cpus=1, input_file=str(sample_yaml_file))
        assert mock_validate.called
    
    def test_input_parser_sets_default_values(self, sample_yaml_file):
        """Test that Input_Parser sets appropriate default values."""
        parser = Input_Parser(cpus=1, input_file=str(sample_yaml_file))
        # Test that some default values are set
        assert hasattr(parser, 'cpus')
        assert parser.cpus == 1


class TestInputParserIntegration:
    """Integration tests for Input_Parser with real YAML files."""
    
    def test_parse_complete_config(self, sample_yaml_config, temp_dir):
        """Test parsing a complete configuration file."""
        yaml_file = temp_dir / "complete_config.yaml"
        with open(yaml_file, 'w') as f:
            yaml.dump(sample_yaml_config, f)
        
        parser = Input_Parser(cpus=2, input_file=str(yaml_file))
        
        # Verify key configuration values are parsed
        assert parser.methodology == 'genetic_algorithm'
        assert parser.population_size == 4
        assert parser.number_of_generations == 2
    
    def test_parse_minimal_config(self, temp_dir):
        """Test parsing a minimal configuration file."""
        minimal_config = {
            'optimization': {
                'methodology': 'genetic_algorithm',
                'population_size': 2,
                'number_of_generations': 1
            }
        }
        
        yaml_file = temp_dir / "minimal_config.yaml"
        with open(yaml_file, 'w') as f:
            yaml.dump(minimal_config, f)
        
        parser = Input_Parser(cpus=1, input_file=str(yaml_file))
        assert parser.methodology == 'genetic_algorithm'
        assert parser.population_size == 2
        assert parser.number_of_generations == 1


@pytest.mark.parametrize("methodology,expected", [
    ("genetic_algorithm", "genetic_algorithm"),
    ("simulated_annealing", "simulated_annealing"),
    ("bayesian_optimization", "bayesian_optimization"),
])
def test_methodology_validation(methodology, expected, temp_dir):
    """Test validation of different optimization methodologies."""
    config = {
        'optimization': {
            'methodology': methodology,
            'population_size': 2,
            'number_of_generations': 1
        }
    }
    
    yaml_file = temp_dir / "methodology_test.yaml"
    with open(yaml_file, 'w') as f:
        yaml.dump(config, f)
    
    parser = Input_Parser(cpus=1, input_file=str(yaml_file))
    assert parser.methodology == expected
