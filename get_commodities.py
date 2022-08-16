import requests

base_currency = 'EUR'

symbol = 'ETH'
endpoint = 'latest'
access_key='l2ka4vji2o0veiovrla2s37u05hqz1rtvi6p81yljqr54mt0qpq15638i6c2'

resp = requests.get(
    'https://commodities-api.com/api/' + endpoint + '?access_key='+ access_key + '&base='+base_currency+'&symbols='+symbol
)

data=resp.json()

print(1/data['data']['rates'][symbol])

