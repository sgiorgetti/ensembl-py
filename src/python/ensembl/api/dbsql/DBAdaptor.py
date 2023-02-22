__all__ = ['ArgumentError']

from ensembl.database.dbconnection import DBConnection

# class DBAdaptor():
#     """Contains all the region related functions over seq_region and coord_system ORM

#     Attributes:
#         dbc: DBConnection
#     """
#     def __init__(self, dbconnection: DBConnection = None) -> None:
#         if not dbconnection:
#             raise Exception(f'Cannot instantiate {self.__class__.__name__}: missing dbconnection')
#         self._dbc = dbconnection

#     def __repr__(self) -> str:
#         """Returns a string representation of this object."""
#         return f'{self.__class__.__name__}({self._dbc._engine})'

#     def get_dbconnection(self):
#         return self._dbc

   
class ArgumentError(ValueError):
    """Exception raised for errors in function calls."""
    def __init__(self, *args: object) -> None:
        super().__init__(*args)