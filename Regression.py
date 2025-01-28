import cupy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LinearRegression
import joblib
# Load data
data_df = pd.read_csv('/home/group5/Data.csv')
label_df = pd.read_csv('/home/group5/Label.csv')

# Transpose data and labels to match samples
X = data_df.values
y = label_df.values

# Split data into training and testing sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Initialize linear regression model
model = LinearRegression()

#Fit the line to the training data
model.fit(X_train, y_train)

import matplotlib.pyplot as plt

# Predictions on the test set
y_pred = model.predict(X_test)

# Plot actual vs. predicted values
plt.figure(figsize=(10, 6))
plt.scatter(y_test, y_pred, color='blue', label='Actual vs. Predicted')
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'k--', lw=2, label='Ideal Prediction')
plt.xlabel('Actual')
plt.ylabel('Predicted')
plt.title('Actual vs. Predicted Values')
plt.legend()
plt.show()
plt.savefig("Regression.png")

mean_test_error = np.mean( (y_test - model.predict(X_test))**2 )

print ('Test MSE: ', mean_test_error)

# Save the trained model to a file
joblib.dump(model, 'Regression.pkl')
