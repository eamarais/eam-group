#######################################################################################################################################
# An example of making "class", "object", "__init__", "methods" in Python
# By defining a class of objects, you can save infomation as different attributes of this object
# so that later you can access the data fields immediately

# Example 1
# Here we define a class of objects, named "User". Then we can save some information of interest as the attributes.

class User:
    def __init__(self,full_name,birthday,language):
    ''' create an object named "User", create some fields (attributes) for it'''
        self.name = full_name
        name_pieces = full_name.split(" ")
        self.first_name = name_pieces[0]
        self.last_name = name_pieces[-1]
        self.birthday = birthday
        self.favourite_language = language

# Now you can build some functions to process data fields associated with this object
# The results can be saved as new attributes of this object
# Or you can "call" the function to access the corresponding results

import datetime

    def age(self):
        '''build a function to calculate the age of this user'''
        today = datetime.date(2020,8,17)
        yyyy = int(self.birthday[0:4])
        mm = int(self.birthday[4:6])
        dd = int(self.birthday[6:8])
        dob = datetime.date(yyyy,mm,dd)
        age_in_days = (today - dob).days
        age_in_years = age_in_days/365
        return age_in_years 

# An example input
user = User("Michael Jordan", "19910101","Python")

# Check the results
print(user.name)
print(user.first_name)
print(user.last_name)
print(user.birthday)
print(user.favourite_language)

# you need to call this function to make it work (since you did not save "age" as an attibute)
print(user.age()) 


# Example 2
# Here we define a class of objects, called "Rectangle". Then we perform some calculations based on its data fields (attributes).
# All functions rely on what has been already provided. For a nested function (e.g. "calculate cost"), it requires the inner functions (e.g. "get_area") to be recognized.

class Rectangle:
    def __init__(self, length, breadth, unit_cost=0):
        self.length = length
        self.breadth = breadth
        self.unit_cost = unit_cost
   
    def get_perimeter(self):
        return 2 * (self.length + self.breadth)
   
    def get_area(self):
        return self.length * self.breadth
   
    def calculate_cost(self):
        area = self.get_area()
        return area * self.unit_cost
# breadth = 120 cm, length = 160 cm, 1 cm^2 = Rs 2000
r = Rectangle(160, 120, 2000)
print("Area of Rectangle: %s cm^2" % (r.get_area()))
print("Cost of rectangular field: Rs. %s " %(r.calculate_cost()))

# End
#######################################################################################################################################
