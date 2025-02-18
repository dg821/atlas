import pandas as pd
import os
from typing import Union, List, Optional, Dict

class ExcelReader:
    """A class to handle reading and processing Excel and CSV files."""
    
    def __init__(self, file_path: str):
        """
        Initialize the ExcelReader with a file path.
        
        Args:
            file_path (str): Path to the Excel or CSV file
        """
        self.file_path = file_path
        self.dataframes: Dict[str, pd.DataFrame] = {}
        self._validate_file()
        self.is_csv = self.file_path.endswith('.csv')
        self.sheet_names = self._get_sheet_names()
        
    def _validate_file(self) -> None:
        """Validate the file exists and has correct extension."""
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"File not found at: {self.file_path}")
            
        if not self.file_path.endswith(('.xlsx', '.xls', '.xlsm', '.csv')):
            raise ValueError("File must be an Excel file (.xlsx, .xls, .xlsm) or CSV file (.csv)")
        
    def _clean_dataframe(self, df: pd.DataFrame) -> pd.DataFrame:
        """Clean the dataframe by removing empty rows/columns and handling missing values."""
        if df.empty:
            return df
            
        # Remove empty rows and columns
        df = df.dropna(how='all')
        df = df.dropna(axis=1, how='all')
        
        return df
    
    def _get_sheet_names(self) -> List[str]:
        """Get all sheet names from the Excel file or return ['Sheet1'] for CSV."""
        if self.is_csv:
            return ['Sheet1']
        try:
            excel_file = pd.ExcelFile(self.file_path)
            return excel_file.sheet_names
        except Exception as e:
            raise Exception(f"Error getting sheet names: {str(e)}")
    
    def read_sheet(self, sheet_name: Union[str, int] = 0, 
                  skiprows: int = 0, 
                  usecols: Optional[List[str]] = None) -> pd.DataFrame:
        """
        Read a specific sheet from the Excel file or the CSV file.
        
        Args:
            sheet_name: Sheet name or index (0-based) - ignored for CSV files
            skiprows: Number of rows to skip at the start
            usecols: List of columns to use (A, B, C, etc.)
            
        Returns:
            pandas DataFrame containing the data
        """
        try:
            if self.is_csv:
                if isinstance(sheet_name, int) and sheet_name > 0:
                    raise ValueError("CSV files only have one sheet")
                df = pd.read_csv(
                    self.file_path,
                    na_values=['NA', 'N/A', ''],
                    keep_default_na=True,
                    skiprows=skiprows,
                    usecols=usecols
                )
            else:
                df = pd.read_excel(
                    self.file_path,
                    sheet_name=sheet_name,
                    engine='openpyxl',
                    na_values=['NA', 'N/A', ''],
                    keep_default_na=True,
                    skiprows=skiprows,
                    usecols=usecols
                )
            
            # Clean the data
            df = self._clean_dataframe(df)
            
            # Store the dataframe
            sheet_key = 'Sheet1' if self.is_csv else (
                sheet_name if isinstance(sheet_name, str) 
                else self.sheet_names[sheet_name]
            )
            self.dataframes[sheet_key] = df
            
            return df
            
        except Exception as e:
            raise Exception(f"Error reading {'CSV' if self.is_csv else 'sheet'}: {str(e)}")
    
    def read_all_sheets(self) -> Dict[str, pd.DataFrame]:
        """Read all sheets from the Excel file or the single CSV sheet."""
        try:
            if self.is_csv:
                self.read_sheet()
            else:
                for sheet_name in self.sheet_names:
                    self.read_sheet(sheet_name)
            return self.dataframes
        except Exception as e:
            raise Exception(f"Error reading all sheets: {str(e)}")

    # Rest of the methods remain the same