import openmeteo_requests
import datetime


class IncreaseSpeed:
  def __init__(self, current_speed: int, max_speed: int):
    self.current_speed = current_speed
    self.max_speed = max_speed
    self.step=10

  def __iter__(self):
    return self
  
  def __next__(self):
    if self.current_speed <= (self.max_speed - self.step):
        self.current_speed += self.step
        return self.current_speed
    raise StopIteration

    
class DecreaseSpeed:
  def __init__(self, current_speed: int):
    self.current_speed = current_speed
    self.step = 10

  def __iter__(self):
    return self
  
  def __next__(self):
    if self.current_speed >= self.step:
      self.current_speed -= self.step
      return self.current_speed
    raise StopIteration
      

class Car:
  cars_on_the_road = 0

  def __init__(self, max_speed: int, current_speed=0, state=True):
    self.max_speed = max_speed
    self.current_speed = current_speed
    self.state = state
    Car.cars_on_the_road += 1
    self.step = 10

  def accelerate(self, upper_border=None):
    if self.state == False:
      print("The car cannot accelerate when it is being parked. Changing the state now")
      self.state = True

    if upper_border == "None":
      if self.current_speed == self.max_speed:
        print("The speed cannot be increased, current speed is: ", self.current_speed)
      elif self.current_speed >= (self.max_speed - self.step):
        self.current_speed == self.max_speed
        print("The speed is increased to maximum, current speed is: ", self.current_speed)
      else:
        self.current_speed += self.step
        print("Current speed is:", self.current_speed)

    elif upper_border != "None" and upper_border <= self.max_speed:
      for i in IncreaseSpeed(self.current_speed, upper_border):
        self.current_speed = i
        print("Current speed is: ", self.current_speed)

      
  def brake(self, lower_border=None):
    if lower_border == "None":
      if self.current_speed < self.step:
        self.current_speed = 0
        print("Current speed is: ", self.current_speed)
      else:
        self.current_speed -= self.step
        print("Current speed is: ", self.current_speed)

    elif lower_border !="None" and lower_border >= 0:
      for i in DecreaseSpeed(self.current_speed, lower_border):
        self.current_speed = i
        print("Current speed is: ", self.current_speed)


  def parking(self):
    if self.state:
      print("The car is currently on the road, going to parking now")
      self.state = False
    else:
      print("The car is already parked")

  @classmethod
  def total_cars(cls):
    print("Total number of cars on the road: ", cls.cars_on_the_road)

  @staticmethod
  def show_weather():
    openmeteo = openmeteo_requests.Client()
    url = "https://api.open-meteo.com/v1/forecast"
    params = {
    "latitude": 59.9386, # for St.Petersburg
    "longitude": 30.3141, # for St.Petersburg
    "current": ["temperature_2m", "apparent_temperature", "rain", "wind_speed_10m"],
    "wind_speed_unit": "ms",
    "timezone": "Europe/Moscow"
    }

    response = openmeteo.weather_api(url, params=params)[0]

    current = response.Current()
    current_temperature_2m = current.Variables(0).Value()
    current_apparent_temperature = current.Variables(1).Value()
    current_rain = current.Variables(2).Value()
    current_wind_speed_10m = current.Variables(3).Value()

    print(f"Current time: {datetime.fromtimestamp(current.Time()+response.UtcOffsetSeconds())} {response.TimezoneAbbreviation().decode()}")
    print(f"Current temperature: {round(current_temperature_2m, 0)} C")
    print(f"Current apparent_temperature: {round(current_apparent_temperature, 0)} C")
    print(f"Current rain: {current_rain} mm")
    print(f"Current wind_speed: {round(current_wind_speed_10m, 1)} m/s")