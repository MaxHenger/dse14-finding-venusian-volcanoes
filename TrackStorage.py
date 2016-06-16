# -*- coding: utf-8 -*-
"""
Created on Mon May 30 13:32:21 2016

@author: MaxHenger
"""

import struct
import numpy as np

def _isArray_(var):
	return isinstance(var, list) or isinstance(var, np.ndarray)

def _isString_(string):
	return isinstance(string, str)

class DataStorage:
	class DataVariable:
		ArrayType = 0
		ScalarType = 1
		StringType = 2

		def __init__(self, name):
			self._name_ = name
			self._values_ = None
			self._axes_ = None
			self._bytes_ = None

		def __toIntConverter__(self, value):
			res = [int(0)] * len(value)

			for i in range(0, len(value)):
				res[i] = int(value[i])

			return res

		def getName(self):
			return self._name_

		def getValues(self):
			return self._values_

		def getNumAxes(self):
			return len(self._axes_)

		def getAxis(self, index):
			return self._axes_[index]

		def setValues(self, values, axes=None):
			self._values_ = values

			# Create header for this specific variable
			stream = bytearray()

			if _isArray_(values):
				values = np.asarray(values)
				numAxes = 0

				if axes != None:
					numAxes = len(axes)

				stream.extend(struct.pack('<BII' + (len(values.shape) * 'I'),
					self.ArrayType, numAxes, len(values.shape), *values.shape))

				if axes != None:
					for iAxis in range(0, len(axes)):
						stream.extend(struct.pack('<I' + str(len(axes[iAxis])) + 'd',
							len(axes[iAxis]), *axes[iAxis]))

				total = values.shape[0]

				for i in range(1, len(values.shape)):
					total *= values.shape[i]

				stream.extend(struct.pack('<' + str(total) + 'd', *values.flatten()))
			elif _isString_(values):
				stringBytes = values.encode('utf-8')
				stream.extend(struct.pack('<BI' + str(len(values)) + 's',
					self.StringType, len(values), stringBytes))
			else:
				stream.extend(struct.pack('<Bd', self.ScalarType, values))

			self._bytes_ = stream

		def setBytes(self, stream):
			# Clear internal data
			self._values_ = []
			self._axes_ = []
			self._bytes_ = []

			# Load data from stream
			self._bytes_ = stream

			offset = struct.calcsize('B')
			dataType = struct.unpack_from('B', stream, 0)[0]

			if dataType == self.ArrayType:
				# Array type, first read header data
				numAxes, shapeSize = struct.unpack_from('<II', stream, offset)
				offset += struct.calcsize('<II')
				shape = struct.unpack_from('<' + str(shapeSize) + 'I', stream, offset)
				offset += struct.calcsize('<' + str(shapeSize) + 'I')

				# Unpack axis data
				for iAxis in range(0, numAxes):
					axisLen = struct.unpack_from('<I', stream, offset)[0]
					offset += struct.calcsize('<I')
					axisData = struct.unpack_from('<' + str(axisLen) + 'd', stream, offset)
					offset += struct.calcsize('<' + str(axisLen) + 'd')

					self._axes_.append(axisData)

				total = shape[0]

				for i in range(1, len(shape)):
					total *= shape[i]

				self._values_ = np.asarray(struct.unpack_from('<' + str(total) + 'd', stream, offset)).reshape(shape)
			elif dataType == self.StringType:
				length = struct.unpack_from('<I', stream, offset)[0]
				offset += struct.calcsize('<I')
				self._values_ = struct.unpack_from('<' + str(length) + 's', stream, offset)[0].decode('utf-8')
				offset += struct.calcsize('<' + str(length) + 's')
			else:
				self._values_ = struct.unpack_from('<d', stream, offset)[0]

	def __init__(self):
		self.variables = []

	def addVariable(self, name, variable, axes=None):
		dataVariable = self.DataVariable(name)
		dataVariable.setValues(variable, axes)
		self.variables.append(dataVariable)

	def getNumVariables(self):
		return len(self.variables)

	def getVariable(self, variable):
		if _isString_(variable):
			for i in range(0, len(self.variables)):
				if self.variables[i]._name_ == variable:
					return self.variables[i]

			return None

		return self.variables[variable]

	def save(self, filename):
		fh = open(filename, 'wb')

		# Write variable header
		stream = bytearray()
		stream.extend(struct.pack('<I', len(self.variables)))

		for i in range(0, len(self.variables)):
			stream.extend(struct.pack('<II' + str(len(self.variables[i]._name_)) + 's',
				len(self.variables[i]._bytes_), len(self.variables[i]._name_),
				self.variables[i]._name_.encode('utf-8')))
			stream.extend(self.variables[i]._bytes_)

		fh.write(stream)
		fh.close()

	def load(self, filename):
		fh = open(filename, 'rb')
		self.variables = []

		numVars = struct.unpack('<I', fh.read(struct.calcsize('<I')))[0]

		for i in range(0, numVars):
			size, stringLength = struct.unpack('<II', fh.read(struct.calcsize('<II')))
			variableName = struct.unpack('<' + str(stringLength) + 's', fh.read(stringLength))[0].decode('utf-8')
			newVariable = self.DataVariable(variableName)
			newVariable.setBytes(fh.read(size))
			self.variables.append(newVariable)

	def info(self):
		print(len(self.variables), 'variables contained:')

		for i in range(0, len(self.variables)):
			print(' >', i + 1, ':', self.variables[i].getName())

			if _isArray_(self.variables[i]._values_):
				print('     type: array')
				print('     shape:', self.variables[i]._values_.shape)
			elif _isString_(self.variables[i]._values_):
				print('     type: string')
				print('     length:', len(self.variables[i]._values_))
			else:
				print('     type: scalar')

			if self.variables[i]._axes_ != None:
				print('     number of axes:', len(self.variables[i]._axes_))

				for j in range(0, len(self.variables[i]._axes_)):
					print('      >', j + 1, ': shape:', self.variables[i]._axes_[j].shape)

			else:
				print('     number of axes: 0')

#ts = DataStorage();
#ts.load('stitched_62000.0at3.5to38000.0at38000.0.dat')
