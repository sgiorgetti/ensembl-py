# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Easy initialisation functionality to set an event logging system."""

__all__ = ["LogLevel", "init_logging"]

import logging
from pathlib import Path
from typing import Optional, Union


LogLevel = Union[int, str]


def init_logging(
    log_level: LogLevel = "WARNING",
    log_file: Optional[Path] = None,
    log_file_level: LogLevel = "DEBUG",
    msg_format: str = "%(asctime)s\t%(levelname)s\t%(message)s",
    date_format: str = r"%Y-%m-%d_%H:%M:%S"
) -> None:
    """Initialises the logging system.

    By default, all the log messages corresponding to `log_level` (and) above will be printed in the
    standard error. If `log_file` is provided, all messages of `log_file_level` level (and above) will
    be written into the provided file.

    Args:
        log_level: Minimum logging level for the standard error.
        log_file: Logging file where to write debug (and above) logging messages.
        log_file_level: Minimum logging level for the logging file.
        msg_format: A format string for the logged output as a whole. More information:
            https://docs.python.org/3/library/logging.html#logrecord-attributes
        date_format: A format string for the date/time portion of the logged output. More information:
            https://docs.python.org/3/library/logging.html#logging.Formatter.formatTime

    """
    # Configure the basic logging system, setting the root logger to the minimum log level available
    # to avoid filtering messages in any handler due to "parent delegation". Also close and remove any
    # existing handlers before setting this configuration.
    logging.basicConfig(format=msg_format, datefmt=date_format, level="DEBUG", force=True)
    # Set the correct log level of the new StreamHandler (by default it is set to NOTSET)
    logging.root.handlers[0].setLevel(log_level)
    if log_file:
        # Create the log file handler and add it to the root logger
        formatter = logging.Formatter(msg_format, datefmt=date_format)
        file_handler = logging.FileHandler(log_file)
        file_handler.setLevel(log_file_level)
        file_handler.setFormatter(formatter)
        logging.root.addHandler(file_handler)
