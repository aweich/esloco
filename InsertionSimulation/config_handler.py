# config_handler.py
import configparser
import ast

class ConfigHandler:
    def __init__(self, config_file="sim_config.ini", section="I"):
        self.config = configparser.ConfigParser()
        self.config.read(config_file)
        self.section = section

    def parse_config(self):
        """Parses the configuration file and returns a dictionary of parameters with correct data types."""
        config_dict = {}
        
        # Iterate through each key-value pair in the specified section
        for key, value in self.config[self.section].items():
            config_dict[key] = self._parse_value(value)
            
        return config_dict

    def _parse_value(self, value):
        """
        Converts a string value to an appropriate Python type.
        """
        if value.lower() == 'none':  # Check for None
            return None
        elif value.lower() == 'true':  # Convert booleans
            return True
        elif value.lower() == 'false':
            return False
        try:
            # Try to convert to an int or float
            return int(value) if '.' not in value else float(value)
        except ValueError:
            pass
        try:
            # Try to parse as dictionary
            parsed_dict = ast.literal_eval(value)
            if isinstance(parsed_dict, dict):
                return parsed_dict
			# Try to parse as a list if it looks like a list
            if ',' in value:
                return [self._parse_value(v.strip()) for v in value.split(',')]
        except (SyntaxError, ValueError):
            pass
        
        # If all else fails, return as string
        return value