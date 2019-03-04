import sys,re

NUM_SEARCH_RANGE = 4
INCREMENT_PERCENT = 0.3
SEARCH_PATTERN = r'<<<(.*)>>>'

class MatchType:

    def __init__(self):
        self.pattern=SEARCH_PATTERN

        self.type = "None"

        self._matched = False
        self._int = 0
        self._float = 0.0
        self._string = "0.0"

    def is_matched(self):
        return self._matched 

    def get_value(self):
        if self.type == 'int':
            return self._int
        elif self.type == 'float':
            return self._float
        else:
            return self._string

    def set_value(self,value,set_type):
        if set_type == 'int':
            self._set_int(value)
        elif set_type == 'float':
            self._set_float(value)
        else:
            self._set_string(value)

    def set_matched(self):
        self._matched = True

    def _set_float(self, value):
        self._float = value
        self.type = 'float'

    def _set_int(self, value):
        self._int = value
        self.type = 'int'

    def _set_string(self, value):
        self._string = value
        self.type = 'string'

def get_num_obj(value):

    m = MatchType()

    match = re.search(m.pattern,value)

    if match:
        value = match.group(1)
        m.set_matched()

    try:
        value = int(value)
        m.set_value(value,'int')
        return m

    except ValueError:

        try: 
            value = float(value)
            m.set_value(value, 'float')
            return m

        except:

            m.set_value(value, 'string')
            return m

def get_parameter_bounds(value):
    increment = INCREMENT_PERCENT*abs(value)
    min_bound = value - NUM_SEARCH_RANGE*increment
    max_bound = value + NUM_SEARCH_RANGE*increment
    return increment, max_bound, min_bound

class ParmGen:

    def __init__(self, ffield):

        self.ffield_name = ffield
        self.param_name = open(ffield+'.param','w')

        self.line = 0
        self.num_section = 0
        self.num_type = 0
        self.num_parameter = 0

        self.header = "" # header

        self.param_general = [] # general
        self.param_elem = [] # element-wise parameters
        self.param_2body = [] # 2body parameters (diagonal)
        self.param_2boff = [] # 2body parameters (off-diagonal)
        self.param_3body = [] # 3 body parameters

        self.section_names = ["none", "general", "element-wise", "2body", "2body offdiag", "angle", "torsion", "hbond"]

        self.fin = open(ffield,'r')

        self.header = self.fin.readline()

        for line in self.fin:
            #print(line, end="")

            for d in line.split():

                m = get_num_obj(d)

                if m.is_matched() and m.type == 'int':
                    self.num_section += 1
                    self.get_param(m.get_value(), self.num_section)

    def get_param(self, num_entry, num_section):
        print('-------------------')
        print('Entering get_param()', self.section_names[num_section], num_entry)

        # skip header lines. element-wise & 2body sections have extra lines
        if num_section == 2:
            self.fin.readline()
            self.fin.readline()
            self.fin.readline()
        elif num_section == 3:
            self.fin.readline()

        for entry in range(num_entry):

            extra_lines = 0
            num_data = 0

            data = self.fin.readline().split()

            # first line of each section stores extra informaiton, such as element name and combination index.
            if num_section == 2:  # element-wise,  e.g. "C    1.3742   4.0000 ..." 
                data = data[1:]
                extra_lines = 3
            elif num_section == 3: # 2body,  e.g. "1  1 141.9346 ..."
                data = data[2:]
                extra_lines = 1
            elif num_section == 4: # 2body off diag, e.g. "1  2   0.0464 ..."
                data = data[2:]
            elif num_section == 5: # angle, e.g. "2  1  2  76.7339 ..."
                data = data[3:]
            elif num_section == 6: # torsion, e.g. "0  3  3  0  -0.9667 ..."
                data = data[4:]
            elif num_section == 7: # hbond, e.g. "3  2  3   2.0000 ..."
                data = data[3:]

            for d in data:
                num_data += 1
                m = get_num_obj(d)

                if m.is_matched():
                    pvalues = get_parameter_bounds(m.get_value())
                    num = [num_section, entry+1, num_data]
                    comment = "%s, %3d-section, %3d-entry, %3d-param, %s"%(self.section_names[num_section], num[0],num[1],num[2],d)
                    self.param_name.write("%4d %4d %4d %10.5f %10.5f %10.5f ! %s\n"%
                        (num[0],num[1],num[2], pvalues[0], pvalues[1], pvalues[2], comment))
                    print("%4d %4d %4d %10.5f %10.5f %10.5f ! %s"%(num[0],num[1],num[2], pvalues[0], pvalues[1], pvalues[2], comment))

            # reading extra lines. no need to handle the combination indices. 
            for num_line in range(extra_lines):
                for d in self.fin.readline().split():
                    num_data += 1
                    m = get_num_obj(d)

                    if m.is_matched():
                        pvalues = get_parameter_bounds(m.get_value())
                        num = [num_section, entry+1, num_data]
                        comment = "%s, %3d-section, %3d-entry, %3d-param, %s"%(self.section_names[num_section], num[0],num[1],num[2],d)
                        self.param_name.write("%4d %4d %4d %10.5f %10.5f %10.5f ! %s\n"%
                            (num[0],num[1],num[2], pvalues[0], pvalues[1], pvalues[2], comment))
                        print("%4d %4d %4d %10.5f %10.5f %10.5f ! %s"%(num_section, entry+1, num_data, pvalues[0], pvalues[1], pvalues[2], comment))


if __name__== "__main__":
    ParmGen(sys.argv[1])
