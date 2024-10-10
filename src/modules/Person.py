

'''
pure infection class 
'''

import numpy as np
class Person:

    # vaccination scenarios ; some default number 
    efficacy = -1
    coverage = -1

    def __init__(self):

        self.health_state = 1
        self.days = 0 # days counter  
        self.expected_exposed_time = -1
        self.expected_infected_time = -1
        self.exposed_date = 0

        # vaccination parameter 
        self.update_vac = False
        #self.ini_days_waned = 0
        self.is_vaccinated = False
        self.protection_level = 0
        
        # day 0 refers to the first day of infection 
        self.vac_periods =  []


    @classmethod
    def get_vac_coverage(cls):
        return cls.coverage
    
    @classmethod
    def get_vac_efficacy(cls):
        if (cls.efficacy == -1):
            raise ValueError
        return cls.efficacy
    
    # def update_vaccination_status(self, vac_tuple ):
    #     '''
    #     in the first year, age of A - B, x% coverage
    #     after the outbreak, every year , vaccinate people who reaches A year old 
    #     we only vaccinate one program once, i.e. for those who are not vaccinated, 
        

    #     AgeOfAdd: the age of the person when he/she is added to the population

    #     update_vac:  if the person has been checked  vaccination before 
    #     is_vaccinated: check if the person is vaccinated 

    #     protection_level: protection level after waning 
    #     '''
        
    #     start_age, vac_time, vac_freq = vac_tuple
        
    #     end_age = start_age + vac_time * vac_freq -1
        
    #     VAC_COVERAGE = self.get_vac_coverage()
    #     VAC_EFFICACY = self.get_vac_efficacy()
        
    #     if(VAC_COVERAGE <0):
    #         raise ValueError
    #         # make it proper 
        
    #     if self.gender == 2 and not(self.update_vac):

    #             # if age of add is within A- B, then  80% prob of getting vaccinated or 
    #             # later she reaches A - B, then she goes into vac group 
    #             # update_vac track if the person is vaccinated                

    #             r1 = np.random.rand()
    #             if  start_age <= self.AgeOfAdd <= end_age:
    #                 self.update_vac = True   
    #                 if r1 <= VAC_COVERAGE:
    #                         self.is_vaccinated = True
    #                         self.protection_level = VAC_EFFICACY 
    #                         vac_total_freq =  vac_freq - (self.AgeOfAdd  - start_age) // vac_time
    #                         vac_total_time =  vac_total_freq * vac_time
    #                         self.vac_periods.append( (0,vac_total_freq * vac_time * 365) ) 
                
    #             elif self.AgeOfAdd  <= start_age and self.age >= start_age :
                    
    #                 self.update_vac = True
                    
    #                 if r1 <= VAC_COVERAGE:
    #                         self.is_vaccinated = True
    #                         self.protection_level = VAC_EFFICACY 
    #                         vac_total_freq =  vac_freq 
    #                         vac_total_time =  vac_total_freq * vac_time
    #                         start_year =   start_age - self.AgeOfAdd # + self.yearOfadd - 2024, as yearOfadd is always 2024
    #                         self.vac_periods.append((start_year * 365,(vac_total_freq * vac_time + start_year) * 365))  
                            

    def infection_direct(self,date):
        '''
        forced infection, ignore the vac effect (almost no effect on final result)
        used for initial infection
        '''
        self.health_state = 3
        self.days = 0
        self.exposed_date = date # assume exposed-date is the infection date 
        self.expected_infected_time = np.round(np.random.normal(5.5,0.77) )
    
    def IfDateWithinPeriod(self,date):
        
        Flag = False

        for vac_period in self.vac_periods:
            if date >= vac_period[0] and \
                date < vac_period[1]:
                    
                Flag = True
            else:
                continue
            
        return Flag
  
    def contact_with_infected(self, date):
        
        r2 = np.random.rand()
        
        if not ( self.IfDateWithinPeriod(date) and r2 <= self.protection_level ):
            # take effect condition means if infection within the vaccination period 
            # resist: the one is vaccinated , date within vaccination, and within the protection level
                self.health_state = 2
                self.exposed_date =date
                # python is always [,), days from day 0 
                self.expected_exposed_time = np.round( np.random.gamma(16.9625,1/2.875) )
                self.expected_infected_time = np.round( np.random.normal(5.5,0.77) ) 
                # self.expected_immune_time = np.random.negative_binomial(1, 1/365/2, size=1)[0] # 365*2 #random.randint(50,365*2 -50)

    def make_infection(self):
        self.health_state = 3
        self.days = 0

    def make_recovery(self):
        self.health_state = 4
        self.days = 0 

    # def back_to_S(self):
    #     # since w assume robost immunity,this function is not usable
    #     self.health_state = 1
    #     self.days = 0
    #     self.expected_exposed_time = -1
    #     self.expected_infected_time  = -1 
    def update_by_day(self):
        self.days += 1


    

if __name__ == "__main__":
    # Go check each step of person instance

    sample_person = Person()
    print("Initial state:", sample_person.__dict__)

    sample_person.infection_direct(30)
    print("Infected directly:", sample_person.__dict__)

    sample_person.contact_with_infected(35)
    print("Contact with infected:", sample_person.__dict__)

    sample_person.update_by_day()
    print("After update by day:", sample_person.__dict__)

    sample_person.make_infection()
    print("After make infection:", sample_person.__dict__)




