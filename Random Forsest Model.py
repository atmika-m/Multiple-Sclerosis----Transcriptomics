!pip install pandas scikit-learn numpy
import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score
import matplotlib.pyplot as plt


file_path = ("")

excel_file = pd.ExcelFile(file_path)
print(excel_file.sheet_names)


deg_df = pd.read_excel(
    file_path,
    sheet_name="SIG_DEGS"
)

print(deg_df.head())



deg_filtered = deg_df[
    (deg_df["P.Value"] <= 0.05) &
    (deg_df["Gene Symbol"].notna())
].copy()


deg_filtered["abs_logFC"] = deg_filtered["logFC"].abs()

deg_filtered = (
    deg_filtered
    .sort_values("abs_logFC", ascending=False)
    .drop_duplicates(subset="Gene Symbol", keep="first")
)


top_up = (
    deg_filtered[deg_filtered["logFC"] > 0]
    .sort_values("logFC", ascending=False)
    .head(100)
)


top_down = (
    deg_filtered[deg_filtered["logFC"] < 0]
    .sort_values("logFC", ascending=True)
    .head(100)
)


selected_genes = (
    top_up["Gene Symbol"].tolist() +
    top_down["Gene Symbol"].tolist()
)

print("Selected genes:", len(selected_genes))


print("Upregulated:", len(top_up))
print("Downregulated:", len(top_down))


data = deg_df.dropna()

target = "logFC" 
X = data[['AveExpr', 't', 'P.Value', 'adj.P.Val','B']]
y = data[target]

valid_rows = y.notna()
X = X.loc[valid_rows]
y = y.loc[valid_rows]

X_train, X_test, y_train, y_test = train_test_split(
    X,
    y,
    test_size=0.2,
    random_state=42
)

rf = RandomForestRegressor(n_estimators=100, random_state=42)
rf.fit(X_train, y_train)
y_pred_rf = rf.predict(X_test)
print("Random Forest Regressor Results (Target =", target, ")")
print("MSE:", mean_squared_error(y_test, y_pred_rf))
print("R2 Score:", r2_score(y_test, y_pred_rf))

lr = LinearRegression()
lr.fit(X_train, y_train)
y_pred_lr = lr.predict(X_test)
print("\nLinear Regression Results (Target =", target, ")")
print("MSE:", mean_squared_error(y_test, y_pred_lr))
print("R2 Score:", r2_score(y_test, y_pred_lr))


from sklearn.ensemble import GradientBoostingRegressor

gb = GradientBoostingRegressor(
    n_estimators=100,
    learning_rate=0.1,
    max_depth=3,
    random_state=42
)
gb.fit(X_train, y_train)
y_pred_gb = gb.predict(X_test)
print("\nGradient Boosting Results (Target =", target, ")")
print("MSE:", mean_squared_error(y_test, y_pred_gb))
print("R2 Score:", r2_score(y_test, y_pred_gb))

plt.figure(figsize=(6,4))
plt.scatter(y_test, y_pred_rf, label="Random Forest", color="red")
plt.scatter(y_test, y_pred_lr, label="Linear Regression", color="blue")
plt.scatter(y_test, y_pred_gb, label="Gradient Boosting", color="green")
plt.plot([min(y_test), max(y_test)], [min(y_test), max(y_test)], 'r--')  # y=x line
plt.xlabel("True Values")
plt.ylabel("Predicted Values")
plt.title(f"Regression Model Predictions (Target = {target})")

plt.legend()
plt.show()


from sklearn.linear_model import LinearRegression
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.tree import DecisionTreeRegressor
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor

models = {
    "Linear Regression": LinearRegression(),
    "Random Forest": RandomForestRegressor(n_estimators=100, random_state=42),
    "Decision Tree": DecisionTreeRegressor(random_state=42),
    "Gradient Boosting": GradientBoostingRegressor(random_state=42),
    "Support Vector Regressor": SVR(kernel="rbf"),
    "KNN Regressor": KNeighborsRegressor(n_neighbors=3)
}
results = {}


for name, model in models.items():
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    results[name] = {"MSE": mse, "R2": r2}
    print(f"{name} → MSE: {mse:.4f}, R²: {r2:.4f}")

plt.figure(figsize=(8,5))
r2_scores = [results[m]["R2"] for m in results]
plt.barh(list(results.keys()), r2_scores, color="blue")
plt.xlabel("R² Score")
plt.title(f"Model Comparison (Target = {target})")
plt.show()


# Predict using trained Linear Regression model
data["Predicted_logFC"] = lr.predict(
    data[['AveExpr', 't', 'P.Value', 'adj.P.Val', 'B']]
)

top100_up = (
    data.sort_values("Predicted_logFC", ascending=False)
        .head(100)
)

print(top100_up[['Gene Symbol','Predicted_logFC']])


top100_down = (
    data.sort_values("Predicted_logFC", ascending=True)
        .head(100)
)

print(top100_down[['Gene Symbol','Predicted_logFC']])

top100_up.to_csv("Top100_Up.csv", index=False)
top100_down.to_csv("Top100_Down.csv", index=False)



data["Predicted_logFC"] = lr.predict(
    data[['AveExpr','t','P.Value','adj.P.Val','B']]
)

data["Residual"] = data["logFC"] - data["Predicted_logFC"]

data["Abs_Residual"] = data["Residual"].abs()
best_genes = (
    data.sort_values("Abs_Residual")
        .head(100)
)

print(best_genes[["Gene Symbol","logFC","Predicted_logFC","Residual"]])

# Predict using trained RF model
data['RF_score'] = rf.predict(
    data[['AveExpr', 't', 'P.Value', 'adj.P.Val','B']]
)

top100_up_rf = (
    data[data['RF_score'] > 0]
    .sort_values(by='RF_score', ascending=False)
    .head(100)
)

display(top100_up_rf)

top100_down_rf = (
    data[data['RF_score'] < 0]
    .sort_values(by='RF_score', ascending=True)
    .head(100)
)

display(top100_down_rf)

data['RF_score'] = rf.predict(
    data[['AveExpr', 't', 'P.Value', 'adj.P.Val','B']]
)

topRF100 = data.sort_values(
    by='RF_score',
    ascending=False
).head(100)
print(topRF100)


top100_up_rf.to_csv("RF_upregulated.csv", index = False)
top100_down_rf.to_csv("RF_downregulated.csv", index = False)
topRF100.to_csv("RF100.csv", index = False)











