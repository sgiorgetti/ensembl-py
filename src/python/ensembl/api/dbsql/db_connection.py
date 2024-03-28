# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""DB connection module for DuckDB"""

__all__ = ["DBConnection", "URL"]

import contextlib
from urllib.parse import urlparse, ParseResult
from typing import TypeVar
import duckdb

URL = TypeVar('URL', str, ParseResult)

class DBConnection():
    def __init__(self, url: URL) -> None:
        if isinstance(URL) == 'str':
            url = urlparse(url)
        self._url = url
        if url.scheme == 'duckdb':
            if not url.path:
                raise ValueError("Invalid URL for duckdb: please, specify ':memory:', ':default:' or a path to file")
        raise NotImplementedError("Can only manage duckdb at the moment")

    def __repr__(self) -> str:
        """Returns a string representation of this object."""
        return f'{self.__class__.__name__}({self._url!r})'

    @property
    def url(self) -> str:
        """Returns a string representation of the connection URL"""
        return self._url
    
    @property
    def dialect(self) -> str:
        """Returns a string representation of the SQL dialect/DB engine"""
        return self._url.scheme

    @contextlib.contextmanager
    def connection(self):
        """Create a connection to the DB specified in URL"""
        if self._url.scheme == "duckdb":
            p = None if self._url.path in (":memory:", ":default:") else self._url.path
            dbc = duckdb.connect(p)
            try:
                yield dbc
                dbc.commit()
            except:
                dbc.rollback()
                raise
            finally:
                dbc.close()
        raise NotImplementedError(f"Cannot connect to DB specified by URL '{self.url}'")
