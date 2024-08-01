"""
Derived from msmc-tools/generate_multihetsep.py
"""
import os
import io
import gzip

class MaskIterator:
	# ZRL: looks like this only supports one forward pass through the bed file
	def __init__(self, filename, negative=False):
		if filename.endswith(".gz"):
			self.file = io.TextIOWrapper(gzip.open(filename, "r"))
		else:
			self.file = open(filename, "r")
		self.eof = False
		self.lastPos = 1
		self.negative = negative
		self.readLine()

	def readLine(self):
		try:
			line = next(self.file)
			fields = line.strip().split()
			if len(fields) == 2:
				self.start = int(fields[0])
				self.end = int(fields[1])
			else:
				self.start = int(fields[1]) + 1
				self.end = int(fields[2])
		except StopIteration:
			self.eof = True
	
	def getVal(self, pos):
		assert pos >= self.lastPos
		self.lastPos = pos
		while pos > self.end and not self.eof:
			self.readLine()
		if pos >= self.start and pos <= self.end:
			return True if not self.negative else False
		else:
			return False if not self.negative else True

class MergedMask:
	def __init__(self, mask_iterators):
		self.maskIterators = mask_iterators

	def getVal(self, pos):
		return all((m.getVal(pos) for m in self.maskIterators))