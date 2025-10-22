---
title: "Analyzing Ocean Wave Buoy Data with Python"
description: "Learn how to download, process, and analyze real ocean wave buoy data using Python. Includes practical examples with NDBC buoy data."
pubDate: 2024-08-05
category: "data-analysis"
tags: ["Python", "buoy data", "NDBC", "time series", "spectral analysis"]
featured: false
---

## Introduction

Ocean wave buoys provide crucial in-situ measurements of wave conditions. In this tutorial, we'll learn how to work with real buoy data from NOAA's National Data Buoy Center (NDBC).

## What Do Wave Buoys Measure?

Modern wave buoys record:
- Significant wave height ($H_s$)
- Peak wave period ($T_p$)
- Mean wave direction
- Directional wave spectrum
- Water temperature
- Air pressure and wind

## Accessing NDBC Data

### Available Data

NDBC provides free access to:
- Real-time data (last 45 days)
- Historical data (archived by year)
- Spectral data
- Continuous winds

### Data Format

NDBC uses a simple text format:

```
#YY  MM DD hh mm  WDIR WSPD  GST  WVHT   DPD   APD MWD  PRES  ATMP  WTMP  DEWP
2024 01 01 00 00  200  12.0  15.0  2.5   8.0   6.5 210 1013.0  15.0  16.0  12.0
```

## Downloading Buoy Data

Let's download data from buoy 46246 (off California):

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

# Buoy ID and year
buoy_id = '46246'
year = 2024

# NDBC URL structure
url = f'https://www.ndbc.noaa.gov/view_text_file.php?filename={buoy_id}h{year}.txt.gz&dir=data/historical/stdmet/'

# Download and parse
df = pd.read_csv(
    url,
    delim_whitespace=True,
    skiprows=[1],  # Skip units row
    na_values=[99.0, 999.0, 9999.0]
)

# Create datetime index
df['datetime'] = pd.to_datetime(df[['#YY', 'MM', 'DD', 'hh', 'mm']])
df = df.set_index('datetime')

print(df.head())
```

## Data Processing

### Handling Missing Data

```python
# Check for missing values
print(df.isnull().sum())

# Fill short gaps (< 3 hours)
df['WVHT'] = df['WVHT'].interpolate(method='time', limit=3)

# Remove long gaps
df = df.dropna(subset=['WVHT'])
```

### Quality Control

```python
def quality_control(df):
    """Apply basic QC to buoy data"""

    # Remove physically impossible values
    df = df[df['WVHT'] >= 0]
    df = df[df['WVHT'] <= 30]  # Max wave height
    df = df[df['DPD'] >= 1]
    df = df[df['DPD'] <= 30]   # Realistic periods

    # Remove statistical outliers (> 5 std)
    for col in ['WVHT', 'DPD']:
        mean = df[col].mean()
        std = df[col].std()
        df = df[np.abs(df[col] - mean) <= 5 * std]

    return df

df = quality_control(df)
```

## Wave Statistics

### Time Series Analysis

```python
# Plot wave height time series
fig, axes = plt.subplots(2, 1, figsize=(12, 8))

# Significant wave height
axes[0].plot(df.index, df['WVHT'], linewidth=0.5)
axes[0].set_ylabel('Wave Height (m)')
axes[0].set_title(f'Buoy {buoy_id} - Wave Conditions')
axes[0].grid(True, alpha=0.3)

# Peak period
axes[1].plot(df.index, df['DPD'], linewidth=0.5, color='orange')
axes[1].set_ylabel('Peak Period (s)')
axes[1].set_xlabel('Date')
axes[1].grid(True, alpha=0.3)

plt.tight_layout()
plt.show()
```

### Statistical Summaries

```python
# Calculate statistics
stats = {
    'Mean Hs': df['WVHT'].mean(),
    'Max Hs': df['WVHT'].max(),
    'Std Hs': df['WVHT'].std(),
    'Mean Tp': df['DPD'].mean(),
    'P90 Hs': df['WVHT'].quantile(0.90),
    'P99 Hs': df['WVHT'].quantile(0.99),
}

print("Wave Statistics:")
for key, value in stats.items():
    print(f"{key}: {value:.2f}")
```

## Advanced Analysis

### Monthly Climatology

```python
# Group by month
monthly = df.groupby(df.index.month).agg({
    'WVHT': ['mean', 'max', 'std'],
    'DPD': 'mean'
})

# Plot monthly statistics
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(monthly.index, monthly['WVHT']['mean'], 'o-', label='Mean')
ax.fill_between(
    monthly.index,
    monthly['WVHT']['mean'] - monthly['WVHT']['std'],
    monthly['WVHT']['mean'] + monthly['WVHT']['std'],
    alpha=0.3
)
ax.set_xlabel('Month')
ax.set_ylabel('Significant Wave Height (m)')
ax.set_title('Monthly Wave Climate')
ax.legend()
ax.grid(True, alpha=0.3)
plt.show()
```

### Storm Detection

```python
def detect_storms(df, threshold=4.0, min_duration_hours=6):
    """Identify storm events based on wave height threshold"""

    # Mark times above threshold
    above_threshold = df['WVHT'] > threshold

    # Find continuous events
    events = []
    in_event = False
    event_start = None

    for idx, above in above_threshold.items():
        if above and not in_event:
            in_event = True
            event_start = idx
        elif not above and in_event:
            event_end = idx
            duration = (event_end - event_start).total_seconds() / 3600

            if duration >= min_duration_hours:
                event_data = df.loc[event_start:event_end]
                events.append({
                    'start': event_start,
                    'end': event_end,
                    'duration_hours': duration,
                    'max_hs': event_data['WVHT'].max(),
                    'mean_hs': event_data['WVHT'].mean()
                })
            in_event = False

    return pd.DataFrame(events)

storms = detect_storms(df, threshold=4.0)
print(f"Detected {len(storms)} storm events")
print(storms)
```

### Wave Rose

```python
from windrose import WindroseAxes

# Create wave rose (waves coming from)
ax = WindroseAxes.from_ax()
ax.bar(df['MWD'], df['WVHT'], normed=True, opening=0.8, edgecolor='white')
ax.set_title('Wave Rose - Direction Distribution')
ax.set_legend()
plt.show()
```

## Spectral Analysis

### Downloading Spectral Data

```python
# Spectral data URL
spec_url = f'https://www.ndbc.noaa.gov/data/realtime2/{buoy_id}.spec'

# Parse spectral data (more complex format)
# Each spectrum has frequency bins and directional information
```

### Computing Wave Parameters

From the spectrum $S(f)$, we can compute:

$$
m_n = \int_0^\infty f^n S(f) \, df
$$

$$
H_s = 4\sqrt{m_0}
$$

$$
T_{m02} = \sqrt{\frac{m_0}{m_2}}
$$

```python
def compute_moments(frequencies, spectrum):
    """Compute spectral moments"""
    df = frequencies[1] - frequencies[0]

    m0 = np.trapz(spectrum, frequencies)
    m1 = np.trapz(frequencies * spectrum, frequencies)
    m2 = np.trapz(frequencies**2 * spectrum, frequencies)

    return m0, m1, m2

# Calculate from spectrum
m0, m1, m2 = compute_moments(freqs, spec)
hs = 4 * np.sqrt(m0)
tm02 = np.sqrt(m0 / m2)
```

## Exporting Results

```python
# Save processed data
df.to_csv(f'buoy_{buoy_id}_{year}_processed.csv')

# Save monthly statistics
monthly.to_csv(f'buoy_{buoy_id}_{year}_monthly.csv')

# Save storm events
storms.to_csv(f'buoy_{buoy_id}_{year}_storms.csv')
```

## Conclusion

NDBC buoy data provides valuable observations for:
- Model validation
- Wave climate studies
- Extreme event analysis
- Engineering design

## Further Reading

- [NDBC Data Description](https://www.ndbc.noaa.gov/docs/ndbc_web_data_guide.pdf)
- [Wave Data Analysis Guide](https://www.ndbc.noaa.gov/wave.shtml)
- Python packages: `pandas`, `xarray`, `windrose`
