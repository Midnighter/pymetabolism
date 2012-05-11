#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
=============================
Reading and Writing KEGG Data
=============================

:Authors:
    Moritz Emanuel Beber
:Date:
    2012-05-11
:Copyright:
    Copyright(c) 2012 Jacobs University of Bremen. All rights reserved.
:File:
    kegg.py
"""


__all__ = ["find_organism"]


import logging
import threading
import re

from .. import miscellaneous as misc

SOAPpy = misc.load_module("SOAPpy", url="http://pywebsvcs.sourceforge.net/")


logger = logging.getLogger(__name__)
logger.addHandler(misc.NullHandler())


class ThreadedKEGGFetcher(threading.Thread):
    """
    docstring
    """
    def __init__(self, wsdl, queue, group=None, target=None, name=None, *args,
            **kw_args):
        """
        docstring
        """
        threading.Thread.__init__(self, group=group, target=target, name=name,
                args=args, kwargs=kw_args)
        # establish connection to DBGET server
        self._serv = SOAPpy.WSDL.Proxy(wsdl)
        self._queue = queue
        self._lock = threading.Lock()

    def run(self):
        """
        docstring
        """
        while True:
            (function, item, output) = self._queue.get()
            try:
                info = eval("self._serv.%s(item)" % function)
            except StandardError:
                logger.debug("psssst:", exc_info=True)
            else:
                if info:
                    self._lock.acquire()
                    output.append(info)
                    self._lock.release()
                else:
                    logger.warn("No information for '%s'.", item)
            finally:
                self._queue.task_done()


def find_organism(self, organism, wsdl="http://soap.genome.jp/KEGG.wsdl",
        browse=10):
    """
    Requires SOAPpy and an active internet connection.
    """
    # establish connection to DBGET server
    serv = SOAPpy.WSDL.Proxy(wsdl)
    # find appropriate organism
    # require user choice here, too many options
    choices = serv.bfind("genome " + organism)
    choices = [choice for choice in choices.split("\n") if choice]
    length = len(choices)
    start = 0
    end = min(length, browse)
    searching = True
    while start < length and searching:
        msg = [""]
        msg.append("Showing organisms %d-%d of %d, please choose an index:"\
                % (start, end - 1, length))
        msg.append("")
        for i in range(start, end):
            msg.append("[%d] %s" % (i, choices[i]))
        if end < length:
            msg.append("")
            msg.append("Type any non-integer to show the next %d organisms."\
                    % min(length - end, browse))
        msg.append("")
        try:
            selection = int(raw_input("\n".join(msg)))
        except ValueError:
            start = end
            end = min(length, end + browse)
        else:
            while True:
                try:
                    choice = choices[selection]
                except IndexError:
                    try:
                        selection = int(raw_input("Chosen index is outside"\
                                " the allowed range, try again:\n"))
                    except ValueError:
                        pass
                else:
                    searching = False
                    break
    logger.info("Please be patient, this will take a few minutes.")
    pattern = re.compile(r"genome:T\d+ (\w+),")
    mobj = pattern.match(choice)
    organism = mobj.group(1)
    return organism

