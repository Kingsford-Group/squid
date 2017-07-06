Nucleotide={'A':'T','C':'G','G':'C','T':'A','R':'Y','Y':'R','S':'W','W':'S','K':'M','M':'K','B':'V','V':'B','D':'H','H':'D', 'N':'N', '.':'.','-':'-'}
def ReverseComp(strs):
	global Nucleotide
	newstrs=""
	for e in strs:
		newstrs+=Nucleotide[e.upper()]
	newstrs=newstrs[::-1]
	return newstrs

class BP_t(object):
	def __init__(self, Chr_, StartPos_, EndPos_, IsLeft_):
		self.Chr=Chr_
		self.StartPos=StartPos_
		self.EndPos=EndPos_
		self.IsLeft=IsLeft_
	def __eq__(self, other):
		if isinstance(other, BP_t):
			return (self.Chr==other.Chr and self.StartPos==other.StartPos and self.EndPos==other.EndPos and self.IsLeft==other.IsLeft)
		return NotImplemented
	def __ne__(self, other):
		result=self.__eq__(other)
		if result is NotImplemented:
			return result
		return not result
	def __lt__(self, other):
		if isinstance(other, BP_t):
			if self.Chr!=other.Chr:
				return self.Chr<other.Chr
			elif self.StartPos!=other.StartPos:
				return self.StartPos<other.StartPos
			elif self.EndPos!=other.EndPos:
				return self.EndPos<other.EndPos
			else:
				return self.IsLeft<other.IsLeft
		return NotImplemented
	def __gt__(self, other):
		if isinstance(other, BP_t):
			if self.Chr!=other.Chr:
				return self.Chr>other.Chr
			elif self.StartPos!=other.StartPos:
				return self.StartPos>other.StartPos
			elif self.EndPos!=other.EndPos:
				return self.EndPos<other.EndPos
			else:
				return self.IsLeft>other.IsLeft
		return NotImplemented
	def __le__(self, other):
		result=self.__gt__(other)
		if result is NotImplemented:
			return result
		return not result
	def __ge__(self, other):
		result=self.__lt__(other)
		if result is NotImplemented:
			return result
		return not result
	def Print(self):
		print("{},{},{},{}".format(self.Chr, self.StartPos, self.EndPos, self.IsLeft))
	def ToString(self):
		return ("{},{},{},{}".format(self.Chr, self.StartPos, self.EndPos, self.IsLeft))
	def ToString(self, RefName):
		return ("{},{},{},{}".format(RefName[self.Chr], self.StartPos, self.EndPos, self.IsLeft))

class SV_t(object):
	def __init__(self, bp1, bp2):
		if bp1<bp2:
			self.BP1=bp1
			self.BP2=bp2
		else:
			self.BP1=bp2
			self.BP2=bp1
	def __eq__(self, other):
		if isinstance(other, SV_t):
			return (self.BP1==other.BP1 and self.BP2==other.BP2)
		return NotImplemented
	def __ne__(self, other):
		result=self.__eq__(other)
		if result is NotImplemented:
			return result
		return not result
	def __lt__(self, other):
		if isinstance(other, SV_t):
			if self.BP1!=other.BP1:
				return self.BP1<other.BP1
			else:
				return self.BP2<other.BP2
		return NotImplemented
	def __gt__(self, other):
		if isinstance(other, SV_t):
			if self.BP1!=other.BP1:
				return self.BP1>other.BP1
			else:
				return self.BP2>other.BP2
		return NotImplemented
	def __le__(self, other):
		result=self.__gt__(other)
		if result is NotImplemented:
			return result
		return not result
	def __ge__(self, other):
		result=self.__lt__(other)
		if result is NotImplemented:
			return result
		return not result
	def Print(self):
		strbp1=self.BP1.ToString()
		strbp2=self.BP2.ToString()
		print(strbp1+"|"+strbp2)
	def ToString(self):
		strbp1=self.BP1.ToString()
		strbp2=self.BP2.ToString()
		return strbp1+"|"+strbp2
	def ToString(self, RefName):
		strbp1=self.BP1.ToString(RefName)
		strbp2=self.BP2.ToString(RefName)
		return strbp1+"|"+strbp2