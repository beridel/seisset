import os
import numpy as np
import datetime as dt
import obspy as obs
from obspy.clients.fdsn import Client
import requests
import ipdb


def get_avail(network, station, location, channels, sampling_rate, dates):
    # stream = requests.get('https://service.iris.edu/fdsnws/availability/1/query?net=UW&sta=DOSE&loc=--&cha=BHZ&nodata=404&format=text', stream=True)
    avail_url = ("https://service.iris.edu/fdsnws/availability/1/query?"
             "net={network}&"
             "sta={station}&"
             "loc={location}&"
             "cha={channel}&"
             "start={starttime}&"
             "end={endtime}&"
             "nodata=404&format=text"
             )

    n_dates = len(dates)
    if not isinstance(channels, list):
        channels = [channels]
    n_channels = len(channels)
    avail = np.zeros((n_dates, n_channels), dtype=np.single)

    start_time = dates[0].strftime('%Y-%m-%dT00:00:00.000000Z')
    stop_time = dates[-1].strftime('%Y-%m-%dT00:00:00.000000Z')

    
    for c in range(n_channels):
        url = avail_url.format(network=network,
                               station=station,
                               location=location,
                               channel=channels[c],
                               starttime=start_time,
                               endtime=stop_time
                               )
        with requests.get(url, stream=False) as avail_stream:
            
            lines = avail_stream.iter_lines(chunk_size=None, decode_unicode=True)
            first_line = next(lines)

            if first_line == 'Error 404: No Data Found':
                break

            # UW DOSE -- BHZ M 40.0 2019-05-03T01:20:59.315000Z 2019-05-03T01:23:10.740000Z
            for line in lines:
                split_line = line.split()
                assert split_line[0] == network
                assert split_line[1] == station
                assert split_line[3] == channels[c]
                
                # who knows what happened, but not worth our time
                if float(split_line[5]) != sampling_rate:
                    continue

                line_start = obs.UTCDateTime(split_line[6])
                line_stop = obs.UTCDateTime(split_line[7])

                line_time = line_start
                line_left = (line_stop - line_start) / 86400 # in days

                # if continuous time is less than 1 hour, skip it; it's not worth our time
                if line_left < 1 / 24:
                    continue

                # check if start and stop time are on the same day
                if line_time.replace(hour=0, minute=0, second=0, microsecond=0) != line_stop.replace(hour=0, minute=0, second=0, microsecond=0):
                    # first day until midnight
                    next_day = (line_time + dt.timedelta(days=1)).replace(hour=0, minute=0, second=0, microsecond=0)

                    day_idx = dates.index(line_time.replace(hour=0, minute=0, second=0, microsecond=0))
                    avail[day_idx, c] += (next_day - line_time) / 86400

                    # now deal with everything after the first day
                    line_time = next_day
                    line_left = (line_stop - line_time) / 86400

                    # can simply add 1 to every day that is complete
                    if line_left > 1:
                        day_idx = dates.index(line_time.replace(hour=0, minute=0, second=0, microsecond=0))
                        avail[day_idx:day_idx + int(line_left), c] += 1
                        
                        line_time = (line_time + dt.timedelta(days=int(line_left))).replace(hour=0, minute=0, second=0, microsecond=0)
                        line_left -= int(line_left)

                # last (partial) day of the line
                day_idx = dates.index(line_time.replace(hour=0, minute=0, second=0, microsecond=0))

                avail[day_idx, c] += (line_stop - line_time) / 86400.

    return avail

all_out_file = 'all_avail.txt'
out_dir = 'avail'
if not os.path.isdir(os.path.join(os.getcwd(), out_dir)):
    os.mkdir(os.path.join(os.getcwd(), out_dir))

networks = 'UW,PB,CC,CN,C8'
min_latitude = 46
max_latitude = 51
min_longitude = -129
max_longitude = -121.5
location = '*'
channels_search = 'BHN,BHE,BHZ,HHN,HHE,HHZ,EHN,EHE,EHZ,EH1,EH2,SHN,SHE,SHZ'
start_time = '2005-01-01T00:00:00'
stop_time = '2023-06-01T00:00:00'
start_before = '2010-01-01T00:00:00'
stop_after = '2020-01-01T00:00:00'
start_date = obs.UTCDateTime(start_time)
stop_date = obs.UTCDateTime(stop_time)
start_before_date = obs.UTCDateTime(start_before)
stop_after_date = obs.UTCDateTime(stop_after)

client = Client('IRIS')
inventory = client.get_stations(network=networks,
                                channel=channels_search,
                                minlatitude=min_latitude,
                                maxlatitude=max_latitude,
                                minlongitude=min_longitude,
                                maxlongitude=max_longitude,
                                starttime=start_date,
                                endtime=stop_date,
                                level='channel')

for network in inventory:
    for station in network:
        earliest_data = stop_date
        latest_data = start_date

        for channel in station:
            if channel.start_date < earliest_data:
                earliest_data = channel.start_date

            if channel.end_date is None:
                latest_data = None
            elif latest_data is not None and channel.end_date > latest_data:
                latest_data = channel.end_date

        if earliest_data > start_before_date:
            inventory = inventory.remove(station=station.code, keep_empty=False)

        if latest_data is not None and latest_data < stop_after_date:
            inventory = inventory.remove(station=station.code, keep_empty=False)

dates = []
date = start_date
while date <= stop_date:
    dates.append(date)
    date += dt.timedelta(days=1)

n_dates = len(dates)

with open(all_out_file, 'w+') as all_out:
    for network in inventory:
        for station in network:
            channel_ids = []
            for channel in station:
                channel_ids.append('{}:{:d}'.format(channel.code, int(channel.sample_rate)))

            channel_ids = set(channel_ids)

            avail_station = np.zeros((n_dates, 3), dtype=np.single) # N E Z

            timer = dt.datetime.now()
            for channel_id in channel_ids:
                channel = channel_id.split(':')[0]
                sampling_rate = float(channel_id.split(':')[1])
                avail = get_avail(network.code, station.code, '*', channel, sampling_rate, dates)

                if channel[-1] == 'N' or channel[-1] == '1':
                    c = 0
                elif channel[-1] == 'E' or channel[-1] == '2':
                    c = 1
                elif channel[-1] == 'Z':
                    c = 2
                avail_station[:, c] += avail.flatten()

                with open(os.path.join(os.getcwd(), out_dir, '{}.{}.{}_avail.txt'.format(network.code, station.code, channel_id)), 'w+') as station_out:
                    for d in range(n_dates):
                        station_out.write('{} {:.2f}\n'.format(dates[d].strftime('%Y-%m-%d'), avail.flatten()[d]))

            if np.sum(avail_station[:, 0:1].flatten()) == 0.:
                out_text = '{}.{} | 1C | shortest: {:.2f} yr | N: {:04.1f} E: {:04.1f} Z: {:04.1f}\n'
            else:
                out_text = '{}.{} | 3C | shortest: {:.2f} yr | N: {:04.1f} E: {:04.1f} Z: {:04.1f}\n'

            channel_uptime = np.sum(avail_station, axis=0)
            nonzero_channels = channel_uptime[channel_uptime !=0]
            if nonzero_channels.size == 0:
                shortest_uptime = 0
            else:
                shortest_uptime = np.min(nonzero_channels) / 365 # in years

            all_out.write(out_text.format(network.code,
                                          station.code,
                                          shortest_uptime,
                                          np.sum(avail_station[:, 0]),
                                          np.sum(avail_station[:, 1]),
                                          np.sum(avail_station[:, 2])))
            print('It took {:.2f}s for {}:{}.'.format((dt.datetime.now() - timer).total_seconds(), network.code, station.code))
